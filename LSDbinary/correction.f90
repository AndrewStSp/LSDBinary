subroutine correction(mode) ! Correct model spectrum corresponding differencies between LSD synthetic model and synthetic spectrum
use observed
use synthetic
use LSD_data
implicit real(8) (a-h,o-z)
real CCF
external CCF_max

 tolerant=1.d-5;  dv=10.d0 ! Code should correct radial velocity by search maximum CCF in region Rv_old +- dv

 select case (mode)
  
  case (1) ! Primary
   Rv1=fmin(Rv1-dv,Rv1+dv,CCF_max,tolerant,CCF,mode)
   allocate (xint(nsynth),yint(nobs))
   xint=wcorr*(1.d0+Rv1/2.997925d5)
   id=map1(xint,corr_prim,nsynth,wave_obs,yint,nobs)
   do i=1,nobs
    dr=yint(i)*(1.d0-model_prim(i))
    model_prim(i)=model_prim(i)-dr
   enddo
   deallocate (xint,yint)

  case (2) ! Secondary
   Rv2=fmin(Rv2-dv,Rv2+dv,CCF_max,tolerant,CCF,mode)
   allocate (xint(nsynth),yint(nobs))
   xint=wcorr*(1.d0+Rv2/2.997925d5)
   id=map1(xint,corr_sec,nsynth,wave_obs,yint,nobs)
   do i=1,nobs
    dr=yint(i)*(1.d0-model_sec(i))
    model_sec(i)=model_sec(i)-dr
   enddo
   deallocate (xint,yint)

 end select   

end

real(8) function CCF_max(v,n)
use synthetic
use LSD_data
implicit real(8) (a-h,o-z)
 
 select case (n) 
  case (1) ! Primary
  allocate (xint(nprim),yint(nprim))
  xint=vLSD_syn1+v
  id=map1(xint,rLSD_syn1,nprim,vLSD_prim,yint,nprim)
  sx=sum(rLSD_prim**2); sy=sum(yint**2); sxy=sum(rLSD_prim*yint)
  if(sx <= 0.or.sy <= 0.or.sxy <= 0)then
   CCF_max=100
  else  
   CCF_max=sqrt(sx)*sqrt(sy)/sxy ! really 1/CCF - fmin search minimum of function, so CCF_max=1/CCF
  endif
  deallocate (xint,yint)
  case (2) ! Secondary
  allocate (xint(nsec),yint(nsec))
  xint=vLSD_syn2+v
  id=map1(xint,rLSD_syn2,nsec,vLSD_sec,yint,nsec)
  sx=sum(rLSD_sec**2); sy=sum(yint**2); sxy=sum(rLSD_sec*yint)
  if(sx <= 0.or.sy <= 0.or.sxy <= 0)then
   CCF_max=100
  else  
   CCF_max=sqrt(sx)*sqrt(sy)/sxy ! really 1/CCF - fmin search minimum of function, so CCF_max=1/CCF
  endif
  deallocate (xint,yint)
 end select

end