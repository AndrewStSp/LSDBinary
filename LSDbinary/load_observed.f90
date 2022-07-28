subroutine load_observed(n) ! Load given observed spectrum
use observed
byte, allocatable :: buf(:) ! arrays for read in FITS format
real(4),allocatable :: fl(:) ! We'll need to have different arrays depending from bitpix !
real(8) w,r,w1,w2,first_pix
character(8) keyword
character(80) string

 w1=wfirst(1); w2=wlast(number_regions)

 select case(if_FITS)

  case (0)  ! ASCII format
  open(10,file=observed_names(n),status='old',iostat=ios)
  if(ios /= 0)then
   write(*,*)trim(observed_names(n))
   write(*,*)'This file not exist'
   stop
  endif
  nobs=0
  do                  ! determine arrays size
   read(10,*,iostat=ios)w
   if(ios /= 0)exit
   if(w < w1)cycle
   if(w > w2)exit
    do i=1,number_regions
    if(w < wfirst(i))cycle
    if(w > wlast(i))cycle
    nobs=nobs+1
   enddo
  enddo  
  allocate (wave_obs(nobs), R_obs(nobs))
  rewind (10)
  nobs=0
  do                  ! determine arrays size
   read(10,*,iostat=ios)w,r
   if(w < w1)cycle
   if(w > w2)exit
   if(ios /= 0)exit
   do i=1,number_regions
    if(w < wfirst(i))cycle
    if(w > wlast(i))cycle
    nobs=nobs+1
	wave_obs(nobs)=w
	R_obs(nobs)=r
   enddo
  enddo
  close (10)

  case (1)  ! FITS format
  open(10,file=observed_names(n),status='old',iostat=ios,access='direct',form='binary',recl=80)
  if(ios /= 0)then
   write(*,*)trim(observed_names(n))
   write(*,*)'This file not exist'
   stop
  endif

 string(1:80)=' '
 nn=0
 do
  nn=nn+1
  read(10,rec=nn)string
  keyword=string(1:8); k=index(string,'/'); if(k == 0)k=81
   select case (keyword)
    case ('BITPIX')
	 read(string(10:k-1),*)bitpix
    case ('CRVAL1')
     read(string(10:k-1),*)w0
    case ('CRPIX1')
     read(string(10:k-1),*)first_pix
    case ('NAXIS1')
	 read(string(10:k-1),*)naxis1
    case ('CDELT1')
	 read(string(10:k-1),*)dwave
    case ('END')
     exit
    end select
 enddo
 length_Header=nn*80  ! Length of Header data
   if(mod(length_Header,2880) == 0)then
      nb=length_Header
   else
      nb=(length_Header/2880+1)*2880
   endif 
  close (10)
  allocate (buf(nb),fl(naxis1))
  open(10,file=observed_names(n),form='binary',convert='BIG_ENDIAN',status='old',iostat=ios)
  read(10)buf,fl
  close (10)

  if(first_pix == 0.d0)first_pix=1.d0
  w0=w0-first_pix*dwave
  w=w0-dwave
  nobs=0
  do i=1,naxis1
   w=w+dwave
   if(w < w1)cycle
   if(w > w2)exit
   do k=1,number_regions
    if(w < wfirst(k))cycle
    if(w > wlast(k))cycle
    nobs=nobs+1
   enddo
  enddo
  
  allocate (wave_obs(nobs), R_obs(nobs))
  w=w0-dwave
  nobs=0
  do i=1,naxis1
   w=w+dwave
   if(w < w1)cycle
   if(w > w2)exit
   do k=1,number_regions
    if(w < wfirst(k))cycle
    if(w > wlast(k))cycle
    nobs=nobs+1
	wave_obs(nobs)=w
	R_obs(nobs)=fl(i)
   enddo
  enddo 
  deallocate (buf,fl)
  
 end select

 R_obs=1.d0-R_obs
 
end     