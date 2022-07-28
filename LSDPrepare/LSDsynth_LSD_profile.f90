subroutine LSD_profile ! Calculate LSD profile
use observed
use LSD_data
use mask
real(8), allocatable :: x(:)
real(8) b
external LSD_fun

if(.not.allocated(model))allocate(model(nobs)) 

allocate (x(nRV))

x = rLSD
call mini(x,R_obs,nRv,nobs,LSD_fun) ! call solver
rLSD = x
deallocate (x)

end

subroutine LSD_fun(x,f,Jacobian,n,nwl,mode) ! calculate model spectrum and Jacobian at each iteration
use LSD_data
use observed
use mask
real(8) x(n), f(nwl), Jacobian(nwl,n), vp, prof, p, dv
real(8), parameter :: c = 2.997925d5  
character(4) mode

 rLSD = x ! LSD profile
 model = 0.d0
 dv = vLSD(2) - vLSD(1)  
 if(mode == 'grad') Jacobian = 0.d0 

 do i = 1, nlines
  do k = 1, nobs 
   vp = (wave_obs(k) - w_mask(i))/w_mask(i)*c 
   if(vp < vLSD(1) .or. vp > vLSD(nRV)) cycle
   id = map1(vLSD,rLSD,nRV,vp,prof,1) 
   model(k) = model(k) + prof*r_mask(i) 
   if(mode == 'grad') then  
    ki = 0
    do l = 1, nRV
	 if(vp <= vLSD(l)) then
	  ki = l; exit
     endif
    enddo
    if(ki == 1) ki = 2; if(ki == 0) ki = nRV
    p = (vp - vLSD(ki-1))/dv 
    Jacobian(k,ki-1) = Jacobian(k,ki-1) + r_mask(i)*(1.d0-p)
    Jacobian(k,ki) = Jacobian(k,ki) + r_mask(i)*p
   endif
  enddo
 enddo

 f = model

end




