subroutine output
use observed
use LSD_data
use mask
real(8), allocatable :: rcorr(:)

 allocate (rcorr(nobs)) 
  
 k = index(observed_spectrum,'.',back=.true.)
 open(300,file = observed_spectrum(1:k)//'lsd')
 open(100,file = observed_spectrum(1:k)//'corr',form='binary')
 rcorr=((1.d0-R_obs)-(1.d0-model))/(1.d0-model)
 do i=1,nRv
   write(300,*)vlsd(i),1-rlsd(i)
 enddo
 close(300)
 write(100)nRv,nobs,vLSD,rLSD,wave_obs,rcorr
 close(100)
 deallocate (rcorr)
 deallocate (wave_obs,R_obs,model)
end
 
