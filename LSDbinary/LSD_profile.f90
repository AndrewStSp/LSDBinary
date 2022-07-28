subroutine LSD_profiles(number) ! Calculate LSD profiles of both components for given observed spectrum
use observed
use LSD_data
real(8), allocatable :: x(:),weights(:),absc(:) ! array of unknowns for solver
real(8) fmin,R2min,R2max
external LSD_fun,R2_fun
real dumf
data niter/2/ ! number of total iterations

tolerant=1.d-5; R2min=0.2d0; R2max=5.d0 ! Diapasone for search square of radia ratio

 allocate (weights(nobs),model_prim(nobs),model_sec(nobs),model(nobs),stat=ios)
 if(ios /= 0)stop 'Memory allocation failed in subroutine LSD_profile'

 weights=1.d0

 model=0.d0 ! Initialize total model spectrum
 call model_spectrum() 

 if(RR2/=0) R2=RR2**2; if(RR2==0) R2=fmin(R2min,R2max,R2_fun,tolerant,dumf,mode)

 nunknowns=nprim+nsec
 allocate (x(nunknowns),absc(nunknowns))

 x(1:nprim)=rLSD_prim; x(nprim+1:nunknowns)=rLSD_sec
 absc(1:nprim)=vLSD_prim; absc(nprim+1:nunknowns)=vLSD_sec

 call mini_reg(x,absc,R_obs,weights,regpar,nunknowns,nobs,LSD_fun) ! call solver
 deallocate (x,absc)

 call output(number)

 deallocate (weights,model_prim,model_sec,model,wave_obs,R_obs)
end

subroutine LSD_fun(x,f,Jacobian,n,nwl,mode) ! calculate model spectrum and Jacobian at each iteration
use LSD_data
use observed
use synthetic
real(8) x(n),f(nwl), Jacobian(nwl,n)
real(8) b,dr,b1,logw
character(4) mode

 rLSD_prim=x(1:nprim); rLSD_sec=x(nprim+1:n) ! LSD profiles of components
 model=0.d0 ! Initialize total model spectrum
 call model_spectrum(Jacobian,mode) ! Calculate model spectra of components and LSD's part of Jacobian
 do i=1,nobs
  logw=log(wave_obs(i))
  b=cp(1)+cp(2)*logw+cp(3)*logw**2 ! Light factor
  b=b*R2
  b1=1.d0/(1.d0+b)
  model(i)=(model_prim(i)+b*model_sec(i))/(1.d0+b) ! Model spectrum
  if(mode == 'grad')then
   do k=1,nprim
    Jacobian(i,k)=Jacobian(i,k)*b1
   enddo
   do k=nprim+1,n
    Jacobian(i,k)=Jacobian(i,k)*b1*b
   enddo
  endif
 enddo 
 f=model
end

real(8) function R2_fun(x,n)
use LSD_data
use observed
use synthetic
real(8) b,s,x,logw
s=0.d0
R2=x
 do i=1,nobs
  logw=log(wave_obs(i))
  b=cp(1)+cp(2)*logw+cp(3)*logw**2 ! Light factor
  b=b*R2
  model(i)=(model_prim(i)+b*model_sec(i))/(1.d0+b) ! Model spectrum
  s=s+(model(i)-R_obs(i))**2
 enddo
 R2_fun=sqrt(s/nobs)
end

