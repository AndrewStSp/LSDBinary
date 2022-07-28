subroutine init_RME ! Initialize variables and read input parameters
use stellar
use synthV_data
character(80) file_name
integer k

 open(1,file='RME.config',status='old',iostat=ios)
 if(ios /= 0)then
  write(999,"(a)")'RME.conf file not exist'; stop
 endif

 read(1,*)ratio ! R2/R1
 ratio=ratio**2
 read(1,*)V_sini !Vsini for primary
 read(1,*)inclin ! inclination angle
 read(1,*)ar1    ! ar=a/R1
 read(1,*)Resolution !Resolved power

 call coordinates

 read(1,'((a))')file_name
 k = index(file_name,'.',back=.true.)
 string_name = trim(adjustl(file_name(1:k-1)))
 call load_synth(1,file_name) ! Load synthetic spectrum of primary
 central_wavelength = (wavelength_last + wavelength_first)/2.d0
 FWHM = central_wavelength/Resolution/(2.d0*sqrt(-log(0.5d0))) ! Full Width at Half Maximum
 FWHM = wavelength_step/FWHM

 read(1,'((a))')file_name
 call load_synth(2,file_name) ! Load synthetic spectrum of secondary
 read(1,*)phase

 close(1,status='delete')

end subroutine init_RME