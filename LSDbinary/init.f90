subroutine init ! Initialize input parameters
!use dfport
use iflport
use observed
use synthetic
use LSD_data
use velocity
character(80) file_mask, file_name
character(5) switch
character(2) ss
character(80) string
character(30) OS_type

 call getenv('OS',os_type)
 slash='\'
 k=index(OS_type,'linux')
 if(k /= 0)slash='/' 

! call getarg(1,file_name) ! Open configuration file with input parameters
 file_name='LSD.config'
 open(1,file=file_name,status='old',iostat=ios)
 if(ios /= 0)stop 'Configuration file not exist'

! Observed spectra info 
 read(1,'((a))')file_mask

 call check_observed_spectra(file_mask)

 read(1,*)number_regions
 allocate(wfirst(number_regions),wlast(number_regions))
 do i=1,number_regions
  read(1,*)wfirst(i),wlast(i)
 enddo

! synthetic spectra info
 read(1,'((a))')synth_path ! path to synthetic spectra data
 k = 0; k = index(synth_path,'#')
 if(k /= 0) synth_path(k:) = ' '

 read(1,'((a))')file_name
 string = trim(synth_path)//trim(file_name)//'.lin'
 call load_mask(1,string); string = '' ! Load mask of primary
 file_prim_corr=trim(synth_path)//trim(file_name)//'.corr'
 
 read(1,'((a))')file_name
 string = trim(synth_path)//trim(file_name)//'.lin'
 call load_mask(2,string); string = '' ! Load mask of secondary
 file_sec_corr=trim(synth_path)//trim(file_name)//'.corr'

 read(1,'((a))')switch
 if(switch == 'table'.or.switch == 'TABLE'.or.switch == 'Table')then
  if_calc=0    ! Read initial guess for radial velocities from file
  read(1,'((a))')file_name
  k = 0; k = index(file_name,'#')
  if(k /= 0) file_name(k:) = ' '
  open(40,file=trim(synth_path)//trim(file_name),status='old',iostat=ios)
  if(ios /= 0)stop 'Cannot open file with radial velocities'
 else
  read(1,*)e ! eccentricity
  if(e == 0.d0)then
   if_calc=1; read(1,*)K1,K2,Vgamma ! circle orbit
  else
   if_calc=2; read(1,*)K1,K2,Vgamma,omega ! elliptical orbit
  endif
  read(1,'((a))')file_name
  open(40,file=trim(synth_path)//trim(file_name),status='old',iostat=ios) 
  if(ios /= 0)stop 'Cannot open file with phases'
 endif
 allocate (f_names(number_spectra),ph(number_spectra),v1(number_spectra),v2(number_spectra),JD(number_spectra))
 do i=1,number_spectra
  f_names(i)= ' '
  read(40,'((a))')string
  k=index(string, ' ')
  f_names(i)(1:k-1)=string(1:k-1)
  if(if_calc == 0)then
   read(string(k:80),*)JD(i),ph(i),v1(i),v2(i)
  else
   read(string(k:80),*)JD(i),ph(i)
  endif
 enddo
 close(40)

 read(1,*)regpar ! Read Tichonov regularization parameters
 read(1,*,iostat = ios) RR2; if(ios /= 0) RR2 = 0
 read(1,*,iostat=ios) ss
 if(ss(1:2)=='E '.or.ss(1:2)=='e '.and.ios==0) then
	read(1,*,iostat = ios) nphases
	if(ios == 0 .and. nphases /= 0) then
		allocate(phases(nphases))
		do i = 1,nphases
			read(1,*)phases(i)
		enddo
	endif
 else
  nphases = 0
 endif
 close(1)
 open(101,file = 'Vrads.dat')
end

subroutine check_observed_spectra(file_mask) ! Check observed spectra format, determine dimensions etc
use observed
!use dflib
use iflport
character(80) file_name,file_mask,string
type (FILE$INFO) info
integer(4) res
integer(kind=int_ptr_kind())handle

 k=index(file_mask,'*')
 if(k == 0)then
  number_spectra=1                   ! single observed file
  allocate(observed_names(1))
  observed_names(1)=trim(file_mask)
 else
  number_spectra=0                   ! multiply observed spectra
  handle = FILE$FIRST
  do
   res = getfileinfoqq(trim(adjustl(file_mask)),info,handle)
   if(handle == FILE$LAST.or.handle == FILE$ERROR) exit
   number_spectra = number_spectra + 1 ! Number of observed spectra spectra
  enddo 
  allocate(observed_names(number_spectra))
  handle = FILE$FIRST; i = 0
  do
   res = getfileinfoqq(trim(adjustl(file_mask)),info,handle)
   if(handle == FILE$LAST.or.handle == FILE$ERROR) exit
   i = i + 1; observed_names(i) = file_mask(1:k-1)//trim(adjustl(info%name)) ! File names of the individual spectra
  enddo
 endif  

! Check observed spectra format
  open(10,file=observed_names(1),status='old',iostat=ios,access='direct',form='binary',recl=80)
  if(ios /= 0)stop 'First observed spectrum file not exist'
  read(10,rec=1)string
  k=index(string,'SIMPLE')
  if(k == 0)then
   if_FITS=0
  else
   if_FITS=1
  endif
  close (10)

end
