program LSDbinary ! 
use observed
use LSD_data
implicit none
integer i

write(*,*)'LSDbinary version 1.01'
 call init
 do i=1,number_spectra
  call load_observed(i)
  call initial_guess(i) ! Get an initial guess for LSD profiles and light factor
  call LSD_profiles(i)
 enddo
 if(allocated(phases)) deallocate(phases)
 stop
end
