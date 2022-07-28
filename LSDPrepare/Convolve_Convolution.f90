subroutine instrumental(R,FWHM,n)  ! Convolve R spectrum by Gaussian with FWHM
 implicit real(8) (a-h,o-z)
 real(8) z(2000),R(n)
 real(8), allocatable :: work(:)
 
  mi=0
  c=FWHM**2
  do i=1,2000
   q=1-i
   d=q*q
   dd=d*c
   if(dd.gt.20.)exit
   mi=mi+1
   z(mi)=exp(-dd)
  enddo
  allocate (work(n))

 do i=1,n
  k1=i-mi
   if(k1 < 1)k1=1
  k2=i+mi
   if(k2 > n)k2=n
  f2=0.d0
  f1=0.d0
   do k=k1,k2
    m=abs(k-i)+1
     if(m > mi)cycle
    f2=f2+z(m)
    f1=f1+R(k)*z(m)
   enddo
   work(i)=f1/f2
  enddo
  R(1:n)=work(1:n)

 deallocate (work)
end

subroutine approximate_rotation(flux,n_sp,Vsini,w0,dw)
 implicit real(8) (a-h,o-z)
 real(8) flux(n_sp)
 real(8), allocatable :: fluxV(:)
 allocate (fluxV(n_sp))
 q=Vsini/2.997925d5
 m1=int(w0*q/(1-q)/dw)+1
 N=n_sp-1
 m2=N-int((w0+dw*n_sp)*q/dw/(1+q))
 if(m2 > n_sp)m2=n_sp
 p2=2.d0/3.1415925d0
 if(m1 > 1)then
  do k=1,m1       ! To avoid end's effects
   fluxV(k)=flux(k)
  enddo
 endif
 if(m2 < n_sp)then
  do k=m2+1,n_sp
   fluxV(k)=flux(k)
  enddo
 endif
  do k=m1+1,m2
   wave=w0+dw*k
   L=int(q*wave/dw)
   z=dw*L/q/wave
   kl1=k-L-1
   if(kl1 < 1)kl1=1
   r1=flux(kL1)
   kl2=k-l
   if(kl2 < 1)kl2=1
   r2=flux(kL2)
   A1=((w0+dw*(kL2))*r1-(w0+dw*(kL1))*r2)/dw
   r3=flux(k+L)
   r4=flux(k+L+1)
   A2=((w0+dw*(k+L+1))*r3-(w0+dw*(k+L))*r4)/dw
   B1=(r2-r1)/dw
   B2=(r4-r3)/dw
   wk=wave
   rck=0.5d0*(a1+a2+wk*(b1+b2))
   sum=fff(k,k+L,-z,flux,w0,dw,q,mode)-fff(k,k-L-1,z,flux,w0,dw,q,mode)
   l2=2*L
    do i=1,l2
     j=k+i-L-1
     ql=1.d0/(q*wk)
     x1=dw*(L-i)*ql
     x2=dw*(L-i+1)*ql
     sum=sum+(fff(k,j,x2,flux,w0,dw,q,mode)-fff(k,j,x1,flux,w0,dw,q,mode))
    enddo
    rck=(rck+p2*sum)
    fluxV(k)=rck
  enddo
 do k=1,n_sp
 flux(k)=fluxV(k)
 enddo
 deallocate (fluxV)
end


real(8) function fff(k,j,x,flux,w0,dw,q,mode)
 implicit real(8) (a-h,o-z)
 real(8) flux(*)	
  if(j < 1)j=1
  wk=w0+dw*k
  wj=w0+dw*j
  if(mode == 1)then
   r1=10**flux(j)
   r2=10**flux(j+1)
  else
   r1=flux(j)
   r2=flux(j+1)
  endif
  aj=((wj+dw)*r1-wj*r2)/dw
  bj=(r2-r1)/dw
  fff=(aj+bj*wk)/2.d0*(x*sqrt(1-x*x)+asin(x))+bj*q*wk*(1-x*x)**1.5d0/3.d0
end 