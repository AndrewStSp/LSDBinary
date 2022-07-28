subroutine disc_integrate(nang,ang, ni,intl,intc, dv, vsi, mac, fluxl,fluxc)
! integrates inten(ni,nang) over ang. Returns flux(ni). 
! Code developed by Christian Stuetz on base of rtint IDL code of J.Valenti

 implicit none

 ! Variables
  integer, parameter :: digi = 8
  real (kind=digi), parameter :: pi = 3.14159265358979323846
  integer          :: ni, nang
  integer          :: ios, n, i, j, k
  real (kind=digi) :: intl(ni,nang), intc(ni,nang), ang(nang)
  real (kind=digi) :: fluxl(ni), fluxc(ni)
  real (kind=digi) :: dv, vsi, mac
  real (kind=digi) :: rmu(nang), sort(nang)
  real (kind=digi) :: r(nang+1), wt(nang)
  integer          :: isort(nang)
  integer (kind=4) :: nm, nmk, nrk
  real (kind=digi) :: sigma, sigr, sigt
  real (kind=digi) :: arg, ar, at
  real (kind=digi) :: norm, v, maxv, r1, r2
  real (kind=digi), allocatable :: rkern(:), ypix(:)
  real (kind=digi), allocatable :: mkern(:), mrkern(:), mtkern(:)

 ! Execute
 
  ! convert input ang to projected radii (rmu) stellar radius is unit
  ! rmu = sin (angle outward normal and line of sight) = sqrt(1-ang)
  do i = 1, nang
    rmu(i) = sqrt(1.0 - ang(i)**2) 
  enddo
  ! sort projected radii and corresponding intensity into ascending order
  ! (i.e. from disk to center of the limb). Equiv. to ang in descending order
  ! make sorted index vector
  do i = 1, nang
   isort(i) = i
  enddo
  do i = 1, nang-1
    do j = i+1, nang
      if ( rmu(i) > rmu(j) ) then
        k = isort(i)
        isort(i) = isort(j)
        isort(j) = k
      endif
    enddo
  enddo
  ! sort rmu()
  do i = 1, nang
    sort(i) = rmu(i)
  enddo
  do i = 1, nang
    rmu(i) = sort(isort(i))
  enddo
  ! sort angles ( ang() )
  do i = 1, nang
    sort(i) = ang(i)
  enddo 
  do i = 1, nang
    ang(i) = sort(isort(i))
  enddo

  ! sort intenitiesni,:)
  do j = 1, ni
    do i = 1, nang
      sort(i) = intl(j,i)
    enddo
    do i = 1, nang
      intl(j,i) = sort(isort(i))
    enddo
    do i = 1, nang
      sort(i) = intc(j,i)
    enddo
    do i = 1, nang
      intc(j,i) = sort(isort(i))
    enddo
  enddo
  ! calculate projected radii for the boundaries of the disk integration
  ! annuli. ( ??? ob dies stimmt )
  r(1) = 0.0d0
  do i = 2, nang
    r(i) = sqrt( 0.5 * ( rmu(i-1)**2 + rmu(i)**2 ) )
  enddo
  r(nang+1) = 1.0d0
  ! integration weights for each disk integration annulus
  do i = 1, nang
    wt(i) = r(i+1)**2 - r(i)**2
  enddo

  ! determine max # of points for rotation velocity kernel
  nrk = 2 * int( vsi*1.0d0/dv ) + 3
  ! determine # of points for the macroturbulence kernel
  sigma = mac / dv / sqrt(2.0)
  nmk = max( 3, min( int(10*sigma), int((ni-3)/2) ) )
  nm = 2 * nmk + 1

  ! allocate arrays
  allocate( rkern(nrk), stat = ios )
    if ( ios /= 0 ) stop "Could not allocate rkern(:) !"
  allocate( mkern(nm), stat = ios )
    if ( ios /= 0 ) stop "Could not allocate mkern(:) !"
  allocate( mrkern(nm), stat = ios )
    if ( ios /= 0 ) stop "Could not allocate mrkern(:) !"
  allocate( mtkern(nm), stat = ios )
    if ( ios /= 0 ) stop "Could not allocate mtkern(:) !"
  i = max( ni+2*nrk , ni+2*nm )
  allocate( ypix(i), stat = ios )
    if ( ios /= 0 ) stop "Could not allocate ypix(:) !"
  ! make sure flux is zero at start
  do i = 1, ni
    fluxl(i) = 0.0
    fluxc(i) = 0.0
  enddo

  ! loop over annuli

  do n = 1, nang
    ! construct the convolution kernel, which describes the distribution
    ! of rotational velocities present in the current annulus. The 
    ! distribution has been derived analytically for annuli of arbitrary
    ! thickness in a rigidly rotating star. There are 2 kernel pieces:
    ! a) radial velocities less than the maximum velocity along the inner
    !    edge of the annulus.
    ! b) velocities bigger than this limit.
    if ( vsi > 0.0d0 ) then
      r1 = vsi * r(n)
      r2 = vsi * r(n+1)
      maxv = r2
      nrk = 2 * int( maxv / dv ) + 3
      norm = 0.0d0
      do i = 1, nrk
        v = abs( dv * ( i - (nrk-1)/2 ) )
        if ( v .lt. r1 ) then                  ! low velocity points
           rkern(i) =   sqrt( r2*r2 - v*v )                   &
          &           - sqrt( r1*r1 - v*v )
           norm = norm + rkern(i)
        elseif ( v .le. r2 ) then             ! high velocity points
           rkern(i) =   sqrt( r2*r2 - v*v )
           norm = norm + rkern(i)
        else
           rkern(i) = 0
        endif
      enddo
      do i = 1, nrk  ! normalize kernel
        rkern(i) = rkern(i) / norm
      enddo
      ! convolve the intensity profile with the rotational velocity
      ! kernel for this annulus. Pad each end of the profile with as
      ! many points as are in the convolution kernel. This reduces
      ! the Fourier fringing.
      if ( nrk .gt. 3 ) then
      ! spectrum
        ! pad
        do i = 1, nrk
          ypix(i) = intl(1,n)
        enddo
        do i = nrk+1, nrk+ni
          ypix(i) = intl(i-nrk,n)
        enddo
        do i = nrk+ni+1, 2*nrk+ni
          ypix(i) = intl(ni,n)
        enddo
        ! convolve
        norm = 0.0
        do i = 1, ni
          norm = 0.0d0
          do j = 1, nrk
            k = nrk + i - int( (nrk-1) / 2 ) + j
            norm = norm + ypix(k) * rkern(j)
          enddo
          intl(i,n) = norm
        enddo
      ! continuum
        ! pad
        do i = 1, nrk
          ypix(i) = intc(1,n)
        enddo
        do i = nrk+1, nrk+ni
          ypix(i) = intc(i-nrk,n)
        enddo
        do i = nrk+ni+1, 2*nrk+ni
          ypix(i) = intc(ni,n)
        enddo
        ! convolve
        norm = 0.0
        do i = 1, ni
          norm = 0.0d0
          do j = 1, nrk
            k = nrk + i - int( (nrk-1) / 2 ) + j
            norm = norm + ypix(k) * rkern(j)
          enddo
          intc(i,n) = norm
        enddo
      endif
    endif
    ! calculate projected sigma for radial and tangential 
    ! velocity distributions.
    sigr = sigma * ang(n)
    sigt = sigma * sqrt( 1.0 - ang(n)**2 )
    ! construct radial macroturb kernel with sigma = sigr
    if ( sigr .gt. 0.0 ) then
      norm = 0.0
      ! build gaussian
      do i = 1, nm
        arg = max( -20.0d0, -0.5*((i-nmk)/sigr)**2 )
        mrkern(i) = exp(arg)
        norm = norm + mrkern(i)
      enddo
      ! normalize
      do i = 1, nm
        mrkern(i) = mrkern(i)/norm
      enddo
    else
      ! build delta function
      do i = 1, nm
        mrkern(i) = 0.0d0
      enddo
      mrkern(nmk) = 1.0d0
    endif
    ! construct tangential macroturb kernel with sigt
    if ( sigt .gt. 0.0 ) then
      norm = 0.0
      ! build gaussian
      do i = 1, nm
        arg = max( -20.0d0, -0.5*((i-nmk)/sigt)**2 )
        mtkern(i) = exp(arg)
        norm = norm + mtkern(i)
      enddo
      ! normalize
      do i = 1, nm
        mtkern(i) = mtkern(i)/norm
      enddo
    else
      ! build delta function
      do i = 1, nm
        mtkern(i) = 0.0d0
      enddo
      mtkern(nmk) = 1.0d0
    endif
    ! sum radial and tangential components weighted by area
    ar = 0.5    ! area covering the radial macroturb
    at = 0.5    ! area covering the tangential macroturb
    do i = 1, nm
      mkern(i) = ar * mrkern(i) + at * mtkern(i)
    enddo
    ! convolve the total flux profiles. Again pad the spectrum
    ! on both ends to prohibit fourier fringing.
                            ! spectrum
    ! pad
    do i = 1, nm
      ypix(i) = intl(1,n)
    enddo
    do i = nm+1, nm+ni
      ypix(i) = intl(i-nm,n)
    enddo
    do i = nm+ni+1, ni+2*nm
      ypix(i) = intl(ni,n)
    enddo
    ! convolve
    do i = 1, ni
      norm = 0.0d0
      do j = 1, nm
        k = nm + i - int( (nm-1) / 2 ) + j
        norm = norm + ypix(k) * mkern(j)
      enddo
      intl(i,n) = norm
    enddo
                            ! continuum
    ! pad
    do i = 1, nm
      ypix(i) = intc(1,n)
    enddo
    do i = nm+1, nm+ni
      ypix(i) = intc(i-nm,n)
    enddo
    do i = nm+ni+1, ni+2*nm
      ypix(i) = intc(ni,n)
    enddo 
    ! convolve
    do i = 1, ni
      norm = 0.0d0
      do j = 1, nm
        k = nm + i - int( (nm-1) / 2 ) + j
        norm = norm + ypix(k) * mkern(j)
      enddo
      intc(i,n) = norm
    enddo

    ! add contribution of the current annulus to the total
    do i = 1, ni
      fluxl(i) = fluxl(i) + intl(i,n)
      fluxc(i) = fluxc(i) + intc(i,n)
    enddo

  enddo

  ! normalize flux
  do i = 1, ni
    fluxl(i) = fluxl(i) / nang
    fluxc(i) = fluxc(i) / nang
  enddo

  ! deallocate velocity kernel vectors
  deallocate(rkern,mkern,mrkern,mtkern,ypix)

end subroutine 

subroutine lininter(ax, ay, bx, by, na, nb, angles)
! Linearily interpolate bx,by on ax,ay for all angles
! ax and bx have to be monotonic increasing

 implicit none

 ! Variables
  integer, parameter :: digi = 8
  integer, parameter :: nlen = 300
  integer          :: i, j, k, ang
  integer          :: na, nb, angles
  real (kind=digi) :: dx, dy
  real (kind=digi) :: ax(na), ay(na,angles)
  real (kind=digi) :: bx(nb), by(nb,angles)

 ! Execute ( derive ay(j:j) )

  do ang = 1, angles

    k = 1 
    do i = 1, na

      do
        if ( k == nb ) exit
        if ( ax(i) < bx(1) ) exit
        if ( ax(i) >= bx(k) .and. ax(i) < bx(k+1) ) exit
        k = k + 1
      enddo

      if ( k == 1 ) then 
        ! if ax < bx(1) extrapolate
        do j = 2, nb
          if ( bx(j) > bx(1) ) exit
        enddo
        dx = bx(j) - bx(1)
        dy = by(j,ang) - by(1,ang)
        ay(i,ang) = by(1,ang) + (ax(i)-bx(1)) * dy/dx      

      elseif ( k == nb ) then
        ! if ax > bx(nb) extrapolate
        do j = nb, 1, -1
          if ( bx(j) < bx(nb) ) exit
        enddo
        dx = bx(nb) - bx(j)
        dy = by(nb,ang) - by(j,ang)
        ay(i,ang) = by(nb,ang) + (ax(i)-bx(nb)) * dy/dx

      else
        ! if ax > bx(1) and ax < bx(nb) interpolate
        j = k + 1
        dx = bx(j) - bx(k)
        dy = by(j,ang) - by(k,ang)
        ay(i,ang) = by(k,ang) + (ax(i)-bx(k)) * dy/dx
      endif

    enddo

  enddo

end subroutine lininter