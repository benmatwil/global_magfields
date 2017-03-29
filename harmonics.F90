module fftw3

  use iso_c_binding

  implicit none

  include 'fftw3.f03'

end module

module harmonics

  use iso_fortran_env, np => real64
  use fftw3

  implicit none

  real(np), parameter :: pi = acos(-1.0_np), fourpi = 4*pi
  
  integer :: lmax
  real(np) :: rmax

  ! grid coordinates
  integer :: nrad, ntheta, nphi
  real(np), dimension(:), allocatable :: rads, thetas, phis
  character(2) :: nsplitstr

  ! coefficients required for calculating harmonics
  real(np), dimension(:,:), allocatable :: coeff1, coeff2, dfact
  
  ! all different global arrays for storing harmonics
  real(np), dimension(:,:,:), allocatable :: plm
  real(np), dimension(:,:,:), allocatable :: qlm, dqlm, qlm_sin
  complex(np), dimension(:,:), allocatable :: blm0
  complex(np), dimension(:,:,:), allocatable :: blm, alm

  ! fft variables
  type(c_ptr) :: plan
  complex(np), dimension(:), allocatable :: fft_in, fft_out

  contains

  ! ################################################################################

  subroutine calc_grids()

    ! Sets up grids in r, theta, phi for final field

    implicit none

    integer, parameter :: nsplit = 16
    integer :: ip, ir, it, isplit

#if fft
    ntheta = (lmax+1)+1
    nphi = 2*(lmax+1)+1
    nrad = 0
    do while (exp(pi*nrad/(lmax+1)/2) < rmax)
      nrad = nrad + 1
    enddo
    ! nrad = ceiling(2*(lmax+1)*log(rmax)/pi)

    allocate(rads(nrad), thetas(ntheta), phis(nphi))

    rads = exp((pi/(lmax+1)/2)*real([(ir, ir=0,nrad-1)], np))
    print*, rads(1), rads(nrad)
    rads = (rmax - 1)*(rads - 1)/(rads(nrad) - 1) + 1 ! re-scale to 1 to rmax
    ! rads = 1.5_np*(real([(ir, ir=0,nrad-1)], np))**3/(nrad-1)**3 + 1 ! cubic grid

    thetas = pi*real([(it, it=0,ntheta-1)], np)/(ntheta-1)
    phis = 2*pi*real([(ip, ip=0,nphi-1)], np)/(nphi-1)
  
#elif analytic

    ! ntheta = (lmax+1)+1
    ! nphi = 2*(lmax+1)+1
    ! nrad = 0
    ! do while (exp(pi*nrad/(lmax+1)/2) < rmax)
    !   nrad = nrad + 1
    ! enddo
    ! ! nrad = ceiling(2*(lmax+1)*log(rmax)/pi)

    ! allocate(rads(nrad), thetas(ntheta), phis(nphi))

    ! rads = exp((pi/(lmax+1)/2)*real([(ir, ir=0,nrad-1)], np))
    ! print*, rads(1), rads(nrad)
    ! rads = (rmax - 1)*(rads - 1)/(rads(nrad) - 1) + 1 ! re-scale to 1 to rmax
    ! ! rads = 1.5_np*(real([(ir, ir=0,nrad-1)], np))**3/(nrad-1)**3 + 1 ! cubic grid

    ! thetas = pi*real([(it, it=0,ntheta-1)], np)/(ntheta-1)
    ! phis = 2*pi*real([(ip, ip=0,nphi-1)], np)/(nphi-1)

    nrad = 0
    do while (exp(pi*nrad/(lmax+1)/2) < rmax)
      nrad = nrad + 1
    enddo
    ntheta = (lmax+1)+1
    nphi = 2*(lmax+1)+1

    nrad = nsplit*nrad-(nsplit-1)
    ntheta = nsplit*ntheta-(nsplit-1)
    nphi = nsplit*nphi-(nsplit-1)

    allocate(rads(nrad), thetas(ntheta), phis(nphi))

    do ir = 0, nrad-1, nsplit
      rads(ir+1) = exp((pi/(lmax+1)/2)*(ir/nsplit))
    enddo
    print*, rads(1), rads(nrad)
    do ir = 2, nrad, nsplit
      do isplit = 0, nsplit-2
        rads(ir+isplit) = ((nsplit-isplit-1)*rads(ir-1) + (isplit+1)*rads(ir+nsplit-1))/nsplit
      enddo
    enddo

    rads = (rmax - 1)*(rads - 1)/(rads(nrad) - 1) + 1

    thetas = pi*real([(it, it=0,ntheta-1)], np)/(ntheta-1)
    phis = 2*pi*real([(ip, ip=0,nphi-1)], np)/(nphi-1)

    nrad = size(rads,1)
    ntheta = size(thetas,1)
    nphi = size(phis,1)

    write(nsplitstr,'(I2.2)') nsplit

#endif

  end

  ! ################################################################################

  subroutine calc_coeffs()

    ! calculates all the coefficients required for both plms and qlms

    implicit none

    integer :: il, im

    allocate(coeff1(0:lmax,0:lmax), coeff2(0:lmax,0:lmax), dfact(0:lmax,0:lmax))
    coeff1 = 0
    coeff2 = 0
    dfact = 0

    !$omp parallel do private(il)
    do im = 0, lmax
      do il = im, lmax
        coeff1(il, im) = sqrt(real(4*il**2 - 1, np))
        coeff2(il, im) = -sqrt(real(2*il + 1, np)*real((il-1)**2 - im**2, np)/(2*il - 3))
        dfact(il, im) = sqrt(real(il**2 - im**2, np))
      enddo
    enddo

  end

  ! ################################################################################

  subroutine calc_plms(xs)

    ! Calculate the Qlms (normalised Plms) for a grid xs
    ! in the range -1 to 1 (cos(theta) = x)

    implicit none

    integer :: il, im, it, ntheta
    real(np) :: xs(:)
    real(np) :: x, y

    if (.not. allocated(coeff1)) call calc_coeffs()

    ntheta = size(xs,1)

    allocate(plm(0:lmax,0:lmax,ntheta))
    plm = 0

    !$omp parallel do private(x, y, im, il)
    do it = 1, ntheta
      x = xs(it)
      y = sqrt(1 - x**2)
      do im = 0, lmax
        if (im == 0) then
          plm(0, im, it) = sqrt(1.0_np/fourpi)
          plm(1, im, it) = sqrt(3.0_np/fourpi)*x
        else
          plm(im, im, it) = -sqrt(real(2*im + 1, np)/(2*im))*y*plm(im-1, im-1, it)
          if (im+1 <= lmax) plm(im+1, im, it) = sqrt(real(2*im + 3, np))*x*plm(im, im, it)
        endif
        do il = im+2, lmax
          if (im < lmax-1) plm(il, im, it) = (x*plm(il-1, im, it)*coeff1(il, im) + &
            plm(il-2, im, it)*coeff2(il, im))/dfact(il, im)
        enddo
      enddo
    enddo

  end

  ! ################################################################################

  subroutine calc_blm_rsun(synmap, nfilter)

    implicit none

    integer :: ia, il, im
    integer, optional :: nfilter
    integer :: ilat, ilon, nlon, nlat
    real(np) :: dt
    real(np), dimension(0:lmax) :: filter
    integer, dimension(0:lmax) :: ls
    complex(np), dimension(:,:), allocatable :: br_blm
    real(np), dimension(:,:), allocatable :: synmap
  
    nlon = size(synmap,1)
    nlat = size(synmap,2)
    
    ia = 10
    do while (exp(-0.25_np*(pi/ia)**2*lmax*(lmax+1)) < 0.1_np)
      ia = ia + 10
    enddo
    ia = ia - 10

    if (present(nfilter)) then
      if (nfilter /= 0) ia = nfilter
    endif
    
    ls = [(il, il=0,lmax)]
    filter = exp(-0.25_np*ls*(ls+1)*(pi/(ia))**2)

    ! not technically calculating the Blms at the correct points... should this be corrected?
    ! may require not using a fft - perhaps quite slow
    allocate(br_blm(0:lmax,nlat))

    allocate(fft_in(nlon), fft_out(nlon))
    plan = fftw_plan_dft_1d(nlon, fft_in,fft_out, fftw_forward,fftw_measure)
    
    !$omp parallel do private(fft_in, fft_out)
    do ilat = 1, nlat
      fft_in = synmap(:,ilat)
      call fftw_execute_dft(plan, fft_in, fft_out)
      br_blm(:,ilat) = fft_out(1:lmax+1)/nlon
    enddo

    call fftw_destroy_plan(plan)
    deallocate(fft_in, fft_out)

    dt = fourpi/nlat

    allocate(blm0(0:lmax, 0:lmax))
    blm0 = 0
    do im = 0, lmax
      do il = im, lmax
        blm0(il, im) = dt*sum(br_blm(im, :)*plm(il, im, :))*filter(il)
      enddo
    enddo

  end subroutine

  ! ################################################################################

  subroutine calc_qlms(theta)

    ! Calculate the Qlms (normalised Plms) for a grid linear in theta
    ! for the range 0 to 2pi

    implicit none

    integer :: il, im, it, ntheta
    real(np) :: theta(:)
    real(np) :: x, y
    real(np) :: c1, c2

    if (.not. allocated(coeff1)) call calc_coeffs()
    
    ! ntheta = 4*(lmax+1)+1
    ntheta = size(theta,1)
    
    allocate(qlm(0:lmax, 0:lmax, ntheta), dqlm(0:lmax, 0:lmax, ntheta))
    qlm = 0
    dqlm = 0
    !$omp parallel do private(x, y, im, il, c1, c2)
    do it = 1, ntheta
      x = cos(theta(it))
      y = sin(theta(it))
      do im = 0, lmax
        if (im == 0) then
          qlm(0, im, it) = sqrt(1.0_np/fourpi)
          dqlm(0, im, it) = 0
          qlm(1, im, it) = sqrt(3.0_np/fourpi)*x
          dqlm(1, im, it) = -sqrt(3.0_np/fourpi)*y
        else
          c1 = -sqrt(real(2*im + 1, np)/(2*im))
          qlm(im, im, it) = c1*y*qlm(im-1, im-1, it)
          dqlm(im, im, it) = c1*(x*qlm(im-1, im-1, it) + y*dqlm(im-1, im-1, it))
          if (im+1 <= lmax) then
            c2 = sqrt(real(2*im + 3, np))
            qlm(im+1, im, it) = c2*x*qlm(im, im, it)
            dqlm(im+1, im, it) = c2*(x*dqlm(im, im, it) - y*qlm(im, im, it))
          endif
        endif
        do il = im+2, lmax
          if (im < lmax-1) then
            qlm(il, im, it) = (x*qlm(il-1, im, it)*coeff1(il, im) + &
              qlm(il-2, im, it)*coeff2(il, im))/dfact(il, im)
            dqlm(il, im, it) = ((x*dqlm(il-1, im, it)-y*qlm(il-1, im, it))*coeff1(il, im) + &
              dqlm(il-2, im, it)*coeff2(il, im))/dfact(il, im)
          endif
        enddo
      enddo
    enddo

    ! add sin dependence for B_phi component
    qlm_sin = qlm
    qlm_sin = 0

    if (ntheta > 1) then
      ! assuming the set of theta end at the poles
      ! !$omp parallel do
      do it = 2, ntheta-1
        qlm_sin(:,1:,it) = qlm(:,1:,it)/sin(theta(it))
      enddo
      qlm_sin(:,1,1) = qlm_sin(:,1,2)
      qlm_sin(:,1,ntheta) = qlm_sin(:,1,ntheta-1)
    else
      qlm_sin(:,1:,1) = qlm(:,1:,1)/sin(theta(1))
    endif

  end

  ! ################################################################################

  subroutine calc_rdep(rads)

    ! extrapolate the blms for each radii

    implicit none

    integer :: ir, il, im, nrad
    real(np), dimension(:) :: rads
    integer, dimension(0:lmax) :: ls
    real(np), dimension(0:lmax) :: r_rsun, r_rmax, div_fact
    real(np), dimension(:,:), allocatable :: rdep_blm, rdep_alm

    nrad = size(rads,1)
    allocate(rdep_blm(0:lmax,nrad), rdep_alm(0:lmax,nrad))
    
    ls = [(il, il=0,lmax)]
    div_fact = ls + 1 + ls*(rmax**(-2*ls - 1))
    
    !$omp parallel do private(r_rsun, r_rmax)
    do ir = 1, nrad
      r_rsun = rads(ir)**(-ls - 2)
      r_rmax = (rads(ir)/rmax)**(2*ls + 1)

      rdep_blm(:,ir) = r_rsun*(ls + 1 + ls*r_rmax)/div_fact
      rdep_alm(:,ir) = r_rsun*(r_rmax - 1)/div_fact
    enddo

    allocate(blm(0:lmax, 0:lmax, nrad))
    ! !$omp parallel do
    do ir = 1, nrad
      blm(:, :, ir) = blm0
    enddo
    !deallocate(blm0)

    alm = blm
    ! !$omp parallel do
    do im = 0, lmax
      blm(:, im, :) = blm(:, im, :)*rdep_blm(:, :)
      alm(:, im, :) = alm(:, im, :)*rdep_alm(:, :)
    enddo

  end

  ! ################################################################################

  subroutine calc_final_field(br, bt, bp)

    implicit none

    integer :: ir, ip, it, il, im, jm
    complex(np) :: mi, expphi
    real(np), dimension(:,:,:), allocatable :: br, bt, bp
    complex(np), dimension(:,:), allocatable :: bbr, bbt, bbp

    allocate(br(nphi, ntheta, nrad), bt(nphi, ntheta, nrad), bp(nphi, ntheta, nrad))

    !$omp parallel workshare
    br = 0
    bt = 0
    bp = 0
    !$omp end parallel workshare

#if analytic

    ! calculating grid analytically
    !$omp parallel do private(ir, it, ip, im, il, mi, jm, expphi)
    do ir = 1, nrad
      print*, ir
      do it = 1, ntheta
        do ip = 1, nphi
          do im = -lmax, lmax
            jm = abs(im)
            mi = cmplx(0.0_np,1.0_np,np)*jm
            expphi = exp(mi*phis(ip))
            if (im >= 0) then
              do il = lmax, abs(im), -1
                br(ip, it, ir) = br(ip, it, ir) + blm(il, jm, ir)*qlm(il, jm, it)*expphi
                bt(ip, it, ir) = bt(ip, it, ir) + alm(il, jm, ir)*dqlm(il, jm, it)*expphi
                bp(ip, it, ir) = bp(ip, it, ir) + mi*alm(il, jm, ir)*qlm_sin(il, jm, it)*expphi
              enddo
            else
              do il = lmax, abs(im), -1
                br(ip, it, ir) = br(ip, it, ir) + conjg(blm(il, jm, ir)*qlm(il, jm, it)*expphi)
                bt(ip, it, ir) = bt(ip, it, ir) + conjg(alm(il, jm, ir)*dqlm(il, jm, it)*expphi)
                bp(ip, it, ir) = bp(ip, it, ir) + conjg(mi*alm(il, jm, ir)*qlm_sin(il, jm, it)*expphi)
              enddo
            endif
          enddo
        enddo
      enddo
    enddo

#elif fft

    allocate(bbr(nphi-1, ntheta), bbt(nphi-1, ntheta), bbp(nphi-1, ntheta))
    
    allocate(fft_in(size(bbr,1)), fft_out(size(bbr,1)))
    plan = fftw_plan_dft_1d(size(fft_in,1), fft_in,fft_out, fftw_backward,fftw_measure)

    !$omp parallel do private(im, mi, il, bbr, bbt, bbp, jm, fft_in, fft_out)
    do ir = 1, nrad
      if (mod(ir, 20) == 0) print*, ir
      bbr = 0
      bbt = 0
      bbp = 0
      do im = lmax, 0, -1
        mi = cmplx(0,im,np)
        jm = im + 1
        do il = lmax, im, -1
          ! sum over l first in preparation for fft
          bbr(jm, :) = bbr(jm, :) + blm(il, im, ir)*qlm(il, im, :)
          bbt(jm, :) = bbt(jm, :) + alm(il, im, ir)*dqlm(il, im, :)
          bbp(jm, :) = bbp(jm, :) + mi*alm(il, im, ir)*qlm_sin(il, im, :)
        enddo
        if (im > 0) then
          ! use the conjugation for negative m rather than direct calculation
          bbr(nphi-im, :) = conjg(bbr(jm, :))
          bbt(nphi-im, :) = conjg(bbt(jm, :))
          bbp(nphi-im, :) = conjg(bbp(jm, :))
        endif
      enddo
      ! fft for each constant theta to get field values
      do it = 1, ntheta
        fft_in = bbr(:,it)
        call fftw_execute_dft(plan, fft_in, fft_out)
        br(1:nphi-1,it,ir) = fft_out
        br(nphi,it,ir) = fft_out(1)

        fft_in = bbt(:,it)
        call fftw_execute_dft(plan, fft_in, fft_out)
        bt(1:nphi-1,it,ir) = fft_out
        bt(nphi,it,ir) = fft_out(1)

        fft_in = bbp(:,it)
        call fftw_execute_dft(plan, fft_in, fft_out)
        bp(1:nphi-1,it,ir) = fft_out
        bp(nphi,it,ir) = fft_out(1)
      enddo
    enddo

    call fftw_destroy_plan(plan)
    deallocate(fft_in, fft_out)

#endif

  end

  ! ################################################################################

  function calc_point(r, theta, phi)

    ! need to have calculated blm(rsun) first
    ! r and theta calculated by sending the correct info to other routines

    implicit none

    real(np) :: r, theta, phi, calc_point(3)
    integer :: il, im, jm
    complex(np) :: mi, expphi
    real(np) :: br, bt, bp

    call calc_qlms([theta])
    call calc_rdep([r])

    br = 0
    bt = 0
    bp = 0

    do im = -lmax, lmax
      jm = abs(im)
      mi = cmplx(0.0_np,1.0_np,np)*jm
      expphi = exp(mi*phi)
      if (im >= 0) then
        do il = lmax, abs(im), -1
          br = br + blm(il, jm, 1)*qlm(il, jm, 1)*expphi
          bt = bt + alm(il, jm, 1)*dqlm(il, jm, 1)*expphi
          bp = bp + mi*alm(il, jm, 1)*qlm_sin(il, jm, 1)*expphi
        enddo
      else
        do il = lmax, abs(im), -1
          br = br + conjg(blm(il, jm, 1)*qlm(il, jm, 1)*expphi)
          bt = bt + conjg(alm(il, jm, 1)*dqlm(il, jm, 1)*expphi)
          bp = bp + conjg(mi*alm(il, jm, 1)*qlm_sin(il, jm, 1)*expphi)
        enddo
      endif
    enddo

    deallocate(qlm, dqlm, qlm_sin, alm, blm)

    calc_point = [br, bt, bp]

  end

end module
