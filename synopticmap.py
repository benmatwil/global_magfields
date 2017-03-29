from __future__ import division
import numpy as np
from astropy.io import fits
import re, glob
import datetime as dt

def get_pole_value(date, pole=None, average=15):
    if pole == 'n':
        headerarg = 'BANDN3'
    elif pole == 's':
        headerarg = 'BANDS3'
    else:
        print("Haven't given either north ('n') or south ('s') to function")

    datelist = [dt.datetime(int(date[:4]), int(date[4:6]), int(date[6:8])) + dt.timedelta(days=x) for x in range(-15, 15)]

    value = np.mean(np.array([fits.open('hmi/fits/mean/hmi.meanpf_720s.{}_000000_TAI.mf_br.fits'.format(day.strftime('%Y%m%d')))[0].header[headerarg] for day in datelist]))
    
    return value

def pole_filler(map, start, fillvalue):
    """
    Function: nan_filler
    Fills the NaNs in a synoptic map with appropriate values

    map : synoptic map array
    start : either 0 or -1
    fill value : float to fill array
    
    Fills the synoptic map map with the value fillvalue starting from
    the latitude with index start. It will stop when it finds a row with no more NaNs
    so you can fill the north and south pole separately for example."""

    row = np.zeros_like(map[:,0])

    if start == 0:
        dir = 1
    elif start == -1:
        dir = -1

    ilat = start
    while np.sum(np.isfinite(map[:,ilat])) < map.shape[0]:
        ilat += dir

    allnan = False
    for k in range(ilat-dir, start-dir, -dir):
        row[:] = map[:,k]
        rowtest = np.isfinite(map[:,k]) == False
        if sum(rowtest) == map.shape[0]:
            allnan = True
        if allnan == True:
            row[:] = fillvalue
        else:
            row[rowtest] = fillvalue
        #print("Filling row",k,"with", fillvalue)
        map[:,k] = row[:]
    
    # fill initial row with the mean value to ensure the pole is smooth
    map[:,start:start+dir:dir] = map[:,start:start+dir:dir].mean()

    return map

def smooth_map(map):
    
    temp = map.copy()

    for ilon in range(-2,map.shape[0]-2):
        for ilat in range(-2,map.shape[1]-2):
            map[ilon,ilat] = ((temp[ilon-1,ilat-1] + temp[ilon+1,ilat-1] + temp[ilon-1,ilat+1] + temp[ilon+1,ilat+1])/8 +
                (temp[ilon+1,ilat] + temp[ilon,ilat+1] + temp[ilon-1,ilat] + temp[ilon,ilat-1])/4 + 3*temp[ilon,ilat]/2)/3

    return map

# Read in Synoptic map and fill/div(B) = 0 appropriately
dir = 'hmi'

for fitsfilename in glob.glob(dir+'/fits/hmi.Mrdailysynframe_720s.*'):

    if fitsfilename.find('daily') > -1:
        # type = 'daily'
        date = re.findall(r"\d{8}", fitsfilename)[0]
    else:
        # type = 'crot'
        date = re.findall(r"\d{4}", fitsfilename)[0]

    fitsfile = fits.open(fitsfilename)
    br1 = fitsfile[0].data.transpose()
    fitsfile.close(fitsfilename)
    br1 = br1.astype(np.float64)
    nlon = br1.shape[0]
    nlat = br1.shape[1]
    lons = np.linspace(0,2*np.pi,nlon+1)
    lats = np.linspace(-1,1,nlat+1)
    lons = (lons[1:] + lons[:-1])/2
    lats = (lats[1:] + lats[:-1])/2
    print('Synoptic map {} of dimensions {} by {}'.format(fitsfilename, nlon, nlat))

    # fill nans
    # The south pole is the 0 index and north pole nlon/-1 index
    # br_north = -30 * 1e-6 * 1e4 # -30 microT to G
    # br_south = 50 * 1e-6 * 1e4 # 50 microT to G
    br_north = get_pole_value(date, pole='n')
    br_south = get_pole_value(date, pole='s')

    print(br_north, br_south)

    br1 = pole_filler(br1, 0, br_south)
    br1 = pole_filler(br1, -1, br_north)
    
    #br1 = smooth_map(br1)

    # make div(B) = 0
    flat = br1.flatten()
    br1 -= np.sum(flat[np.argsort(abs(flat))])/nlon/nlat

    # write to file in fortran order
    with open(dir+'/synmap_' + date + '.dat', 'wb') as savefile:
        np.array([nlon],dtype=np.int32).tofile(savefile)
        np.array([nlat],dtype=np.int32).tofile(savefile)
        for i in range(nlat):
            br1[:,i].tofile(savefile) # fortran order
        lons.tofile(savefile)
        lats.tofile(savefile)
