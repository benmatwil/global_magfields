import numpy as np
import matplotlib.pyplot as plt
plt.ion()

dir = 'data5v/'
filename = dir+'field_20100601_0351.dat'

with open(filename, 'rb') as file:
  nrad = np.fromfile(file, count=1, dtype=np.int32)
  ntheta = np.fromfile(file, count=1, dtype=np.int32)
  nphi = np.fromfile(file, count=1, dtype=np.int32)
  br = np.fromfile(file, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nrad,ntheta,nphi)
  bt = np.fromfile(file, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nrad,ntheta,nphi)
  bp = np.fromfile(file, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nrad,ntheta,nphi)
  rads = np.fromfile(file, count=nrad, dtype=np.float64)
  thetas = np.fromfile(file, count=ntheta, dtype=np.float64)
  phis = np.fromfile(file, count=nphi, dtype=np.float64)

plt.figure()
plt.contourf(phis, thetas, br[0,:,:], 
  np.linspace(-10,10,21), cmap='RdBu_r', extend='both')

plt.figure()
plt.contourf(phis, thetas, bt[0,:,:], 
  np.linspace(-10,10,21), cmap='RdBu_r', extend='both')
plt.figure()
plt.contourf(phis, thetas, bp[0,:,:], 
  np.linspace(-10,10,21), cmap='RdBu_r', extend='both')

filename = 'data/field_20100601_0081.dat'
br, bt, bp, rads, thetas, phis = rd.field(filename)

plt.figure()
plt.contourf(phis, thetas, br[0,:,:],
  np.linspace(-10,10,21), cmap='RdBu_r', extend='both')
plt.ylim([thetas[-1], thetas[0]])
plt.gca().invert_xaxis()

plt.figure()
plt.contourf(phis, thetas, bt[0,:,:], 
  np.linspace(-10,10,21), cmap='RdBu_r', extend='both')
plt.figure()
plt.contourf(phis, thetas, bp[0,:,:], 
  np.linspace(-10,10,21), cmap='RdBu_r', extend='both')

with open('hmi/synmap_20100601.dat', 'rb') as file:
  nlon = np.fromfile(file, count=1, dtype=np.int32)
  nlat = np.fromfile(file, count=1, dtype=np.int32)
  synmap = np.fromfile(file, count=nlon*nlat, dtype=np.float64).reshape(nlat,nlon)
  lons = np.fromfile(file, count=nlon, dtype=np.float64)
  lats = np.fromfile(file, count=nlat, dtype=np.float64)

plt.figure()
plt.contourf(lons, 180*np.arccos(lats[::-1])/np.pi, synmap, 
  np.linspace(-10,10,21), cmap='RdBu_r', extend='both')

filename = 'data/magfield_0351_20100601.dat'
with open(filename, 'rb') as file:
  nrad = np.fromfile(file, count=1, dtype=np.int32)
  ntheta = np.fromfile(file, count=1, dtype=np.int32)
  nphi = np.fromfile(file, count=1, dtype=np.int32)
  br = np.fromfile(file, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nphi,ntheta,nrad).T
  bt = np.fromfile(file, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nphi,ntheta,nrad).T
  bp = np.fromfile(file, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nphi,ntheta,nrad).T
  rads = np.fromfile(file, count=nrad, dtype=np.float64)
  thetas = np.fromfile(file, count=ntheta, dtype=np.float64)
  phis = np.fromfile(file, count=nphi, dtype=np.float64)

bridl = br[0,:,:]

filename = '/sraid6v1/bmw/global/data/field_20100601.dat'
with open(filename, 'rb') as file:
  nrad = np.fromfile(file, count=1, dtype=np.int32)
  ntheta = np.fromfile(file, count=1, dtype=np.int32)
  nphi = np.fromfile(file, count=1, dtype=np.int32)
  br = np.fromfile(file, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nrad,ntheta,nphi)
  bt = np.fromfile(file, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nrad,ntheta,nphi)
  bp = np.fromfile(file, count=nrad*ntheta*nphi, dtype=np.float64).reshape(nrad,ntheta,nphi)
  rads = np.fromfile(file, count=nrad, dtype=np.float64)
  thetas = np.fromfile(file, count=ntheta, dtype=np.float64)
  phis = np.fromfile(file, count=nphi, dtype=np.float64)

brfortran = br[0,:,:]

absbr = abs(bridl[::2,:] - brfortran)
for i in xrange(brdiff.shape[0]):
  for j in xrange(brdiff.shape[1]):
    a = np.min([abs(bridl[i,j]),abs(brfortran[i,j])])
    brdiff[i,j] = abs(bridl[i,j] - brfortran[i,j])/a



br, bt, bp, rads, thetas, phis = rd.field('data/field_20100601_0081.dat')

for it, t in enumerate(thetas):
  br[:,it,:] = br[:,it,:]*np.sin(t)

ints = np.trapz(np.trapz(br, phis), thetas)
plt.plot(rads,ints)




polar60 = np.zeros(len(files), dtype=np.float64)
polar70 = polar60.copy()
polar60smooth = polar60.copy()
polar70smooth = polar60.copy()

for i, file in enumerate(files):
  f = fits.open(file)
  polar60[i] = f[0].header['CAPN2']
  polar70[i] = f[0].header['BANDN3']
  f.close

for i in xrange(5,len(files)-5):
  polar60smooth[i] = np.mean(polar60[i-5:i+5])
  polar70smooth[i] = np.mean(polar70[i-5:i+5])

pixels = np.zeros(len(files), dtype=np.int32)
for i, file in enumerate(files):
  f = fits.open(file)
  pixels[i] = f[0].header['NUMN3']
  f.close