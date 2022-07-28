subroutine LSDSynth(vsin) ! Program that computes LSD profile from a single star synthetic spectrum
use observed
use mask
real :: vsin

Vsini = vsin

call init 
call LSD_profile 
call output

deallocate(wfirst,wlast,w_mask,r_mask,el_mask)

end subroutine LSDSynth