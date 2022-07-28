subroutine init 
use observed
use res
implicit none
integer(4) :: ios, k ,i

 open(1,file='LSD.conf',status='old',iostat=ios)
 if(ios /= 0)then
  write(999,"(a)")'LSD.conf file not exist'; stop
 endif

 read(1,'((a))') observed_spectrum; k = 0; k = index(observed_spectrum,'!',back=.true.); if(k /= 0) observed_spectrum(k:) = ' '

 read(1,*) number_regions
 allocate(wfirst(number_regions),wlast(number_regions))
 do i = 1,number_regions
  read(1,*) wfirst(i), wlast(i)
 enddo

 call load_observed 

 read(1,'((a))') line_list; k = 0; k = index(line_list,'!',back=.true.); if(k /= 0) line_list(k:) = ' '

 call load_mask 

 close(1,status='delete')

 call initial_guess 

end
