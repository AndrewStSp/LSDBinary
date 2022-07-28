subroutine load_synth(mode,file_name)  ! Load synthetic spectrum, mode =1: local, mode=2: spot
use synthV_data
use coordinates_data
implicit real(8) (a-h,o-z)
character(80) file_name
real(8) Inu_c(7),Inu_l(7),Ii,mu
real(8) b(ndeg+1),sspoly(ndeg+1),stats(10)   ! arrays for approximation by mu
real(8) :: mu_i(7) = (/0.9636d0,0.8864d0,0.8018d0,0.7071d0,0.5976d0,0.4629d0,0.2673d0/)
real(8), allocatable :: b1(:),b2(:),b3(:),b4(:),b5(:),b6(:),b7(:),b8(:),x(:)
parameter (nmu=7)
save


 select case(mode) 
  case (1)     ! read synthetic spectrum of primary and expand it in mu and wavelengths space
  open(10,file=file_name,form='binary',access='direct',recl=144,status='old',iostat=ios)
  if(ios /= 0)stop 'File with synthetic spectrum of primary not exist'

  read(10,rec=1)wavelength_first,wavelength_last,wavelength_step,number_of_wavelengths
  nw=number_of_wavelengths*10
  allocate(I_c(nw,20),I_l(nw,20),b1(number_of_wavelengths),b2(number_of_wavelengths),b3(number_of_wavelengths),stat=i1)
  allocate (b4(number_of_wavelengths),b5(number_of_wavelengths),b6(number_of_wavelengths),x(number_of_wavelengths),stat=i2)
  allocate (b7(number_of_wavelengths),b8(number_of_wavelengths),f1(number_of_wavelengths),f1c(number_of_wavelengths),stat=i3)
  if(i1+i2+i3 /= 0) stop 'Memory allocation failed, too large synthetic spectrum'

   wv=wavelength_first-wavelength_step
   do i=1,number_of_wavelengths
    wv=wv+wavelength_step
	x(i)=wv
    read(10,rec=i+1)flux_c,flux,Inu_c,Inu_l
    call DRCURV(nmu,mu_i,Inu_l,ndeg,b,sspoly,stats)! approximate total intensity
	b1(i)=b(1); b2(i)=b(2); b3(i)=b(3); b4(i)=b(4)
    call DRCURV(nmu,mu_i,Inu_c,ndeg,b,sspoly,stats)! approximate continuum intensity
	b5(i)=b(1); b6(i)=b(2); b7(i)=b(3); b8(i)=b(4)
   enddo
   dw=wavelength_step/10.d0   ! Interpolate local intensities to extended wavelengths and mu sets
   wv=wavelength_first-dw
   do i=1,nw
    wv=wv+dw
	id=map1(x,b1,number_of_wavelengths,wv,b(1),1)
	id=map1(x,b2,number_of_wavelengths,wv,b(2),1)
	id=map1(x,b3,number_of_wavelengths,wv,b(3),1)
	id=map1(x,b4,number_of_wavelengths,wv,b(4),1)
	mu=-dmu
	do mp=1,20
	 mu=mu+dmu
	 Ii=0.d0
     do ip=1,ndeg+1   ! Restore intensity
      Ii=Ii*mu+b(ndeg+2-ip)
     enddo
     I_l(i,mp)=II
    enddo 
   enddo
   wv=wavelength_first-dw
   do i=1,nw
    wv=wv+dw
	id=map1(x,b5,number_of_wavelengths,wv,b(1),1)
	id=map1(x,b6,number_of_wavelengths,wv,b(2),1)
	id=map1(x,b7,number_of_wavelengths,wv,b(3),1)
	id=map1(x,b8,number_of_wavelengths,wv,b(4),1)
	mu=-dmu
	do mp=1,20
	 mu=mu+dmu
	 Ii=0.d0
     do ip=1,ndeg+1   ! Restore intensity
      Ii=Ii*mu+b(ndeg+2-ip)
     enddo
     I_c(i,mp)=II
    enddo 
   enddo

  deallocate (b1,b2,b3,b4,b5,b6,b7,b8,x)
  close (10)

  case (2) ! Read synthetic spectrum of secondary: spectrum was convolved with rotation and instrumewntal profile before
  allocate (f2(number_of_wavelengths),f2c(number_of_wavelengths),stat=i1)
  if(i1 /= 0) stop 'Memory allocation failed, too large synthetic spectrum'

  open(10,file=file_name,status='old',iostat=ios)
  if(ios /= 0)stop 'File with synthetic spectrum of primary not exist'
   do i=1,number_of_wavelengths
    read(10,*)d,d,f2(i),f2c(i)
   enddo
 end select

 close (10)

end