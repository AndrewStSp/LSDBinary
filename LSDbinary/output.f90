subroutine output(n) ! outputs for n observed spectrum
use observed
use LSD_data
use velocity

 k=index(observed_names(n),'.',back=.true.)
 open(100,file=observed_names(n)(1:k)//'lsd1')
 do i=1,nprim  ! debug only: call output
  write(100,*)vLSD_prim(i),1.d0-rLSD_prim(i)
 enddo
 close (100)
 open(100,file=observed_names(n)(1:k)//'lsd2')
 do i=1,nsec
  write(100,*)vLSD_sec(i),1.d0-rLSD_sec(i)
 enddo
 close (100)
 open(100,file=observed_names(n)(1:k)//'prof')
 do i=1,nobs
  write(100,'(f11.4,3f10.3)')wave_obs(i),1.d0-model(i),1.d0-model_prim(i),1.d0-model_sec(i)
 enddo
 close (100)
 write(101,'(f16.6,4f10.3)')BJD,phase,Rv1,Rv2,sqrt(R2)

end
 
