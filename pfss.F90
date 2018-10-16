program pfss_merge

#if _OPENMP
  use omp_lib
#endif

  use iso_fortran_env, only : np => real64
  use harmonics

  implicit none

  ! reading in from the command line
  character(50) :: arg
  integer :: iarg, iinput
  integer :: ic = 0, lc = 0, nfilter = 0, ac = 0, dc = 0, ec = 0
  integer :: nproc
  real(np) :: modelrad = 0
  
  ! timing
  real(np) :: start, finish

  ! reading in synoptic map
  character(50) :: synfilename
  integer :: nlon, nlat
  real(np), allocatable :: synmap(:, :)
  real(np), dimension(:), allocatable :: lons, lats

  ! final field
  real(np), dimension(:, :, :), allocatable :: br, bt, bp, brw, btw, bpw
  character(100) :: outfname, outputdir
  character(4) :: lmaxstr
  character(7) :: alphastr, dstr
  character(:), allocatable :: typestr, fname_sub, calc_type, extrafname, endfname

  ! loop parameters
  integer :: ilat, ilon, ir, ip, it
  
  ! ################################################################################
  
  call cpu_time(start)
  outputdir = 'data'
  rmax = 2.5_np
  zeroalpha = .false.
  nproc = 0

  if (command_argument_count() == 0) then
    stop 'No command line parameters given'
  else
    do iarg = 1, command_argument_count()
      call get_command_argument(iarg, arg)
      if (arg(1:2) == '-i') then
        ! read in synoptic map data file name
        call get_command_argument(iarg+1, arg)
        synfilename = trim(arg)
        ic = 1
      elseif (arg(1:2) == '-o') then
        ! read in output directory
        call get_command_argument(iarg+1, arg)
        outputdir = trim(arg)
      elseif (arg(1:2) == '-l') then
        ! read in lmax
        call get_command_argument(iarg+1, arg)
        read(arg, *) lmax
        lc = 1
      elseif (arg(1:2) == '-f') then
        ! read in nfilter
        call get_command_argument(iarg+1, arg)
        read(arg, *) nfilter
      elseif (arg(1:2) == '-r' .or. arg(1:6) == '--rmax') then
        ! read in model radius R_max
        call get_command_argument(iarg+1, arg)
        read(arg, *) rmax
      elseif (arg(1:2) == '-a' .or. arg(1:7) == '--alpha') then
        ! read in alpha
        call get_command_argument(iarg+1, arg)
        read(arg, *) alpha
        ac = 1
        if (trim(arg) == '0' .or. trim(arg) == '0.' .or. trim(arg) == '0.0' .or. alpha == 0.0_np) then
          zeroalpha = .true.
        endif
      elseif (arg(1:2) == '-d') then
        ! read in d
        call get_command_argument(iarg+1, arg)
        read(arg, *) d
        dc = 1
      elseif (arg(1:2) == '-e') then
        ! read extra filename ending
        call get_command_argument(iarg+1, arg)
        extrafname = trim(arg)
        ec = 1
      elseif (arg(1:4) == '--np') then
        ! read number of processors to use
        call get_command_argument(iarg+1, arg)
        read(arg, *) nproc
      endif
    enddo
    if (ic == 0 .or. lc == 0) then
      stop "Haven't given a data input filename (-i) and lmax (-l)"
    endif
#if mhs
    if (ac == 0 .or. dc == 0) then
      stop "Haven't given a data input alpha (-a or --alpha) and d (-d)"
    endif
#endif
    if (mod(lmax, 2) == 0) stop 'lmax is not odd'
  endif

#if _OPENMP
  if (nproc > 0) call omp_set_num_threads(nproc)
  print*, 'Using', omp_get_max_threads(), 'processors'
#endif

  print*, 'lmax is', lmax
  print*, 'rmax is', rmax
#if mhs
  print*, 'alpha is', alpha
  print*, 'd is', d
#endif

  call calc_grids()

  print*, 'r from', rads(1), 'to', rads(nrad)
  print*, 'theta from', thetas(1), 'to', thetas(ntheta)
  print*, 'phi from', phis(1), 'to', phis(nphi)
  
  print*, 'Number of points', nrad, ntheta, nphi

  ! reading in transformed synoptic map
  open(unit=1, file=synfilename, access='stream', status='old')
    read(1) nlon, nlat
    allocate(synmap(nlon, nlat))
    allocate(lons(nlon), lats(nlat))
    read(1) synmap, lons, lats
  close(1)

  ! calculate coefficients needed for calc_plms and calc_qlms
  print*, 'Calculating coefficients'
  call calc_coeffs()

  print*, 'Calculating Plms'
  call calc_plms(lats)

  print*, 'Calculating Blms for the solar surface'
  call calc_blm_rsun(synmap, nfilter)

  deallocate(plm)

  print*, 'Calculating Qlms'
  call calc_qlms(thetas)

  print*, 'Calculating r dependence for the Blms'
  call calc_rdep(rads)

  deallocate(coeff1, coeff2, dfact, blm0)

  call cpu_time(finish)
  print*, ceiling(finish-start), 'seconds since start'
  call cpu_time(start)
  
  ! now to calculate final field
  print*, 'Calculating final field'
  call calc_final_field(br, bt, bp)

  deallocate(qlm, qlm_sin, dqlm, blm, alm)

  call cpu_time(finish)
  print*, real(finish-start, 4), 'seconds to make grid'

  ! writing final field to file
  write(lmaxstr,'(I4.4)') lmax
  fname_sub = synfilename(index(synfilename, '/',.true.)+len('synmap_')+1:index(synfilename, '.dat')-1)
  if (ec == 1) then
    endfname = '-' // extrafname // '.dat'
  else
    endfname = '.dat'
  endif
#if analytic
  calc_type = 'anal' // nsplitstr
#elif fft
  calc_type = 'fft'
#endif
#if pfss
  outfname = trim(outputdir) // '/newpfss_field_' // &
    fname_sub // '_' // lmaxstr // '_' // calc_type // endfname
#elif mhs
  write(alphastr, '(F7.3)') alpha
  write(dstr, '(F7.3)') d

  if (alpha < 0) then
    if (len(trim(adjustl(alphastr))) == 6) alphastr = '-0'//alphastr(3:7)
  else
    if (len(trim(adjustl(alphastr))) == 5) alphastr = '00'//alphastr(3:7)
    if (len(trim(adjustl(alphastr))) == 6) alphastr = '0'//alphastr(2:7)
  endif

  if (d < 0) then
    if (len(trim(adjustl(dstr))) == 6) dstr = '-0'//dstr(3:7)
  else
    if (len(trim(adjustl(dstr))) == 5) dstr = '00'//dstr(3:7)
    if (len(trim(adjustl(dstr))) == 6) dstr = '0'//dstr(2:7)
  endif
  
  typestr = 'finite'

  outfname = trim(outputdir) // '/newmhs_field_' // &
    fname_sub // '_' // lmaxstr // '_' // calc_type // &
    '_alpha_' // alphastr // '_d_' // dstr // '-' // typestr // endfname
#endif
  print*, 'Writing to file ' // outfname

  allocate(brw(nrad, ntheta, nphi), btw(nrad, ntheta, nphi), bpw(nrad, ntheta, nphi))

  print*, 'Transposing arrays for writing'
  do ip = 1, nphi
    do it = 1, ntheta
      do ir = 1, nrad
        brw(ir, it, ip) = br(ip, it, ir)
        btw(ir, it, ip) = bt(ip, it, ir)
        bpw(ir, it, ip) = bp(ip, it, ir)
      enddo
    enddo
  enddo

  open(unit=2, file=trim(outfname), access='stream', status='replace')
    write(2) size(brw, 1), size(brw, 2), size(brw, 3)
    write(2) brw, btw, bpw
    write(2) rads, thetas, phis
  close(2)

end program
