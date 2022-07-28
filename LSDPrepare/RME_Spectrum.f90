subroutine spectrum ! Calculate spectrum of primary and composite spectrum
use coordinates_data
use synthV_data
use stellar
use eclips
implicit real(8) (a-h,o-z)
character(7) fn

 covering = 0
 pi4=1.d0/(4.d0*pi) 

 write(fn,'(f7.4)')phase

 open(20,file=trim(string_name)//'_'//trim(adjustl(fn))//'.dat',iostat=ios)

!open(21,file=trim(fn)//'_composite.dat',iostat=ios)


 rad = 2.d0*pi*phase
 cosi=cos(inclin*2.d0*pi/360.d0)
 x0=ar1*sin(rad); y0=ar1*cos(rad)*cosi

 wave=wavelength_first-wavelength_step; dw=wavelength_step/10.d0; n10=number_of_wavelengths*10
 do k=1,number_of_wavelengths
  wave=wave+wavelength_step
  flux=0.d0; fluxc=0.d0
  do ix=1,nbin
   x=x_gau(ix)
   wv=wave*z(ix) ! Doppler shift
   nw=(wv-wavelength_first)/dw+1; if(nw > n10)nw=n10; if(nw <= 0)nw=1 ! wavelength index in I arrays
   xp=x+x0
   fl=0.d0; fc=0.d0
   do iy=1,nbin
    y=y_gau(ix,iy)
    n_mu=mu_s(ix,iy)/dmu+1; if(n_mu > 20)n_mu=20 ! mu index in I arrays
    yp=y+y0
    if(xp**2+yp**2 > ratio)then
	 fl=fl+I_l(nw,n_mu)*wt_y(ix,iy)
	 fc=fc+I_c(nw,n_mu)*wt_y(ix,iy)
    endif
   enddo
   if(fl == 0.d0)cycle

   flux=flux+fl*wt_gau(ix)
   fluxc=fluxc+fc*wt_gau(ix)
  enddo
  flux=flux*pi4 ! Transform physical flux in Eddington flux (same as in output of synthetic spectrum)
  fluxc=fluxc*pi4
  f1(k)=flux; f1c(k)=fluxc
 enddo
 if(sum(f1) == 0.d0 .and. sum(f1c) == 0.d0) then
    deallocate (I_c,I_l,f1,f1c,f2,f2c)
    close (20,status='delete')
    covering = 1
    return
 endif
 call instrumental(f1,FWHM,number_of_wavelengths)  ! Convolve with instrumental profile (Gaussian with FWHM)
 call instrumental(f1c,FWHM,number_of_wavelengths)

 wave=wavelength_first-wavelength_step
 do k=1,number_of_wavelengths
  wave=wave+wavelength_step
  write(20,'(f12.3,f10.5,1p2e14.3)')wave,f1(k)/f1c(k),f1(k),f1c(k)
 enddo  
 close (20)
 deallocate (I_c,I_l,f1,f1c,f2,f2c)
! call composite ! calculate the composite spectrum
end




subroutine composite !Produce composed spectrum of binary system
use synthV_data
use stellar
implicit real(8) (a-h,o-z)
real(8), allocatable :: wsynt(:),fsum(:),fcsum(:),wint(:),fint(:)
real(8), parameter :: c=2.997925d5

  allocate (wsynt(number_of_wavelengths),fsum(number_of_wavelengths),fcsum(number_of_wavelengths),wint(number_of_wavelengths),fint(number_of_wavelengths),stat=i1)
  if(i1 /= 0) stop 'Memory allocation failed, too large synthetic spectrum'

 wave=wavelength_first-wavelength_step
 do k=1,number_of_wavelengths
  wave=wave+wavelength_step
  wsynt(k)=wave
 enddo
  nsp=number_of_wavelengths
  wint=wsynt*(1.d0+V1/c)
  id=map1(wint,f1,nsp,wsynt,fint,nsp)
  fsum=fint
  id=map1(wint,f1c,nsp,wsynt,fint,nsp)
  fcsum=fint
  wint=wsynt*(1.d0+V2/c)
  id=map1(wint,f2,nsp,wsynt,fint,nsp)
  fsum=fsum+ratio*fint
  id=map1(wint,f2c,nsp,wsynt,fint,nsp)
  fcsum=fcsum+ratio*fint

 do k=1,number_of_wavelengths
  write(21,'(f12.3,f10.5,1p2e14.3)')wsynt(k),fsum(k)/fcsum(k),fsum(k),fcsum(k)
 enddo  
 close (21)

 deallocate (wsynt,fsum,fcsum,wint,fint)

end
