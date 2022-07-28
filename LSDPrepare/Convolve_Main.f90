subroutine convolve ! Convolution synthetic spectra from synthV by rotation,macroturbulence and Gassian profile
 implicit real(8) (a-h,o-z)
 parameter (nang=7)        ! Number of rings on stellar disk
 real(8) Inu_c(7),Inu_l(7) ! Specific intensities: continuum and total
 character(80) file_name        !  files name
 real(8), allocatable :: flux(:),flux_c(:),I_c(:,:),I_l(:,:),tau_eff(:),numbers_eff(:) !Fluxes,intensities,effective depths,correspondend atmosphere layers,work array
 real(8), allocatable :: wlc(:), spxc(:,:),spxl(:,:),x(:),fxl(:),fxc(:)
 real(8) :: mu(7) = (/0.9636d0,0.8864d0,0.8018d0,0.7071d0,0.5976d0,0.4629d0,0.2673d0/)
 
 open(1,file='convolve.conf')  ! Fixed name configuration file
 read(1,'((a))')file_name 
 open(19,file=file_name ,form='binary',access='direct',recl=144,status='old',iostat=ios)
 if(ios /= 0)stop 'Input synthetic spectrum file not exist'
 read(1,*)Vsini         ! Rotation velocity
 read(1,*,iostat=ios)Resolved_power ! Resolution
 if(ios /= 0)Resolved_power=-1
 read(1,'((a))')file_name 
 open(11,file=file_name)
 Vmacro = 0.d0
 read(19,rec=1)wavelength_first,wavelength_last,wavelength_step,number_of_wavelengths,Teff,glog,V_micro 
 allocate (flux(number_of_wavelengths),flux_c(number_of_wavelengths),I_c(number_of_wavelengths,nang),I_l(number_of_wavelengths,nang),&
           tau_eff(number_of_wavelengths),numbers_eff(number_of_wavelengths),wlc(number_of_wavelengths),stat=j1)
 if(j1 /= 0)stop 'Memory allocation 1 failed'

! Convolution with rotation & macroturbulence
   wl=wavelength_first-wavelength_step
   do i=1,number_of_wavelengths
    wl=wl+wavelength_step; wlc(i)=wl
    read(19,rec=i+1)flux_c(i),flux(i),Inu_c,Inu_l,effective_depth_tau,effective_depth_number
	I_l(i,1:nang)=Inu_l(1:nang); I_c(i,1:nang)=Inu_c(1:nang)
    tau_eff(i)=10**effective_depth_tau; numbers_eff(i)=effective_depth_number
   enddo

  if(vsini /= 0.d0)then   ! Rotation & macroturbulence
   resol=1.d6
   nx = int( log( wavelength_last / wavelength_first ) / log( 1.d0 + 1.d0 / resol ) + 1) ! Interpolate into log spaced scale
   resol = 1.d0 / ( ( wavelength_last / wavelength_first ) ** ( 1.d0 / ( nx - 1.d0 ) ) - 1.d0 )
   dv = 2.99792458d5 / resol ! scale to velocity space
   allocate (x(nx),spxc(nx,nang),spxl(nx,nang),fxl(nx),fxc(nx),stat=j2)
   if(j2 /= 0)stop 'Memory allocation 2 failed'
   x(1) = wlc(1)
   do i = 2, nx
    x(i) = x(i-1) * ( 1.d0 + 1.d0 / resol )
   enddo
   call lininter(x, spxl, wlc, I_l, nx, number_of_wavelengths, nang)
   call lininter(x, spxc, wlc, I_c, nx, number_of_wavelengths, nang)
   call disc_integrate(nang,mu,nx,spxl,spxc, dv, Vsini, Vmacro, fxl,fxc)
   id=map1(x,fxl,nx,wlc,flux,number_of_wavelengths)
   id=map1(x,fxc,nx,wlc,flux_c,number_of_wavelengths)
   deallocate (x,spxc,spxl,wlc,fxl,fxc)
   flux=flux/4.d0; flux_c=flux_c/4.d0  !!! look like normalization during disk intehration let additional factor 4!!!
   call approximate_rotation(tau_eff,number_of_wavelengths,Vsini,wavelength_first,wavelength_step) ! Approximate rotation by Lyubimkov approach 
   call approximate_rotation(numbers_eff,number_of_wavelengths,Vsini,wavelength_first,wavelength_step)
  endif

if(Resolved_power > 0)then
 central_wavelength=( wavelength_last + wavelength_first )/2.d0
 FWHM=central_wavelength/Resolved_power/(2.d0*sqrt(-log(0.5d0))) ! Full Width at Half Maximum
 FWHM=wavelength_step/FWHM
 call instrumental(flux,FWHM,number_of_wavelengths)  ! Convolve with instrumental profile (Gaussian with FWHM)
 call instrumental(flux_c,FWHM,number_of_wavelengths)
 call instrumental(tau_eff,FWHM,number_of_wavelengths)
 call instrumental(numbers_eff,FWHM,number_of_wavelengths)
endif
   do i=1,number_of_wavelengths
    write(11,'(f10.4,f10.5,1p2e15.4,0p2f8.3)')wavelength_first+(i-1)*wavelength_step,flux(i)/flux_c(i),flux(i),flux_c(i),log10(tau_eff(i)),numbers_eff(i)
   enddo
 deallocate (I_l,I_c,flux,flux_c,tau_eff,numbers_eff) 
 close(1,status='delete'); close(11); close(19) 
end