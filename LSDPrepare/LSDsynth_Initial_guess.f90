subroutine initial_guess ! Prepare RV grid for the LSD profile
use observed
use LSD_data
use init_prepare
implicit real(8) (a-h,o-z)
real(8) gx(2), gauss

 nRv = npoint_lsd 
 if(Vsini <= 80.d0)then
  call gauss_profile(gx) 
  vmax = sqrt(-gx(2)*gx(2)*log(min_prof/gx(1)))
  if(vmax < 10)nRv = npoint_lsd/3
 else
  central_wavelength = (wfirst(1) + wlast(number_regions))/2.d0
  call Rotational_kernel(central_wavelength,FWHM_vsini) 
  vmax = FWHM_vsini
 endif
 allocate (vLSD(nRv),rLSD(nRv))
 dv = 2.d0*vmax/float(nRv)
 v = -vmax-dv
 do i = 1,nRv
  v = v+dv
  vLSD(i) = v
  if(Vsini <= 80.d0)then
   rLSD(i) = gauss(v,gx)
  else
   rLSD(i) = 0.1d0
  endif
 enddo
 return
end


real(8) function gauss(v,gx)
 implicit none
 real(8) :: gx(2), v
    gauss = gx(1)*exp(-(v/gx(2))**2)
end function gauss


subroutine gauss_profile(gx) ! Calculate gauss LSD profile
use observed
implicit real(8) (a-h,o-z)
real(8) gx(2)
external gauss_fun
 gx(1) = 1.d0
 gx(2) = 3.d0 
 call mini(gx,R_obs,2,nobs,gauss_fun) 
end


subroutine gauss_fun(gx,f,Jacobian,n,nwl,mode) 
use observed
use mask
use init_prepare
implicit real(8) (a-h,o-z)
real(8) gx(n), f(nwl), Jacobian(nwl,n)
real(8), parameter :: c = 2.997925d5 
character(4) mode

 if(gx(1) < 0.d0) gx(1) = abs(gx(1))
 vmax = sqrt(-gx(2)*gx(2)*log(min_prof/gx(1))); vmin = -vmax 
 f = 0.d0  
 if(mode == 'grad') Jacobian = 0.d0 
 do i = 1, nlines 
  do k = 1, nwl 
   vp = (wave_obs(k) - w_mask(i))/w_mask(i)*c 
   if(vp < vmin .or. vp > vmax)cycle
   f(k) = f(k)+gx(1)*exp(-(vp/gx(2))**2)*r_mask(i) 
   if(mode == 'grad') then 
     const2 = vp/(gx(2))**2
	const = exp(-const2*vp)
	const3 = 2.d0*gx(1)*const*const2
	Jacobian(k,1) = Jacobian(k,1)+ const
	Jacobian(k,2) = Jacobian(k,2)+ const3*vp
   endif
  enddo
 enddo
end


subroutine Rotational_kernel(central_wavelength,FWHM_vsini) 
use observed
integer imax(1)
real(8) central_wavelength, y, w, w1, delta, Gmax, Gmin, Ghalf, FWHM_vsini
real(8), allocatable :: lambda(:), G(:)
real(8), parameter :: c = 2.997925d5, pi = 3.14159d0, epsilon = 0.5d0 
delta = central_wavelength*vsini/c
w1 = int(delta) + 1.d0
m = 0; w = -w1 - wavelength_step
do
 w = w + wavelength_step
 if (w > w1) exit
 y = 1.d0 - (w/delta)**2.d0; if(y < 0.d0) cycle
 m = m + 1
enddo
allocate(lambda(m),G(m))
m = 0; w = -w1 - wavelength_step
do
 w = w + wavelength_step
 if (w > w1) exit
 y = 1.d0 - (w/delta)**2.d0; if(y < 0.d0) cycle
 m = m + 1
 lambda(m) = w
 G(m) = (2.d0*(1.d0-epsilon)*dsqrt(y)+pi*epsilon/2.d0*y)/(pi*delta*(1.d0-epsilon/3.d0)) ! kernel
enddo
Gmax = maxval(G); Gmin = minval(G); imax = maxloc(G)
Ghalf = (Gmax - Gmin)*0.5d0
do i = imax(1), 1, -1
 if(G(i) < Ghalf) exit
enddo
i = i + 1; FWHM_vsini = abs(2.d0*lambda(i)/central_wavelength*c) ! FWHM of the rotational kernel in km/s
deallocate(lambda,G)
end
