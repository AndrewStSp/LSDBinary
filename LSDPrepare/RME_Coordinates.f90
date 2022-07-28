subroutine coordinates ! Calculate cartesian coordinates on visible disk
use coordinates_data
use stellar
implicit real(8) (a-h,o-z)
real(8) mu,yt(nbin),wyt(nbin)
  call MakeGaussMatrix(-1.d0,1.d0,x_gau,wt_gau,nbin) ! Calculate x&weigths for Gauss quadratures

  do ix=1,nbin
   x=x_gau(ix)
   Vrad=x*V_sini ! radial velocity for given x 
   z(ix)=1.d0-Vrad/2.997925d5 ! Doppler factor. Sing - from usage in disc integrate: with red shift we need use blue shifted wavelength point in fixed wavelengths
   y1=sqrt(1.d0-x*x)
   call MakeGaussMatrix(-y1,y1,yt,wyt,nbin)
   y_gau(ix,1:nbin)=yt; wt_y(ix,1:nbin)=wyt
    do iy=1,nbin
     y=yt(iy)
	 p=1.d0-(x*x+y*y)
	 mu=0.d0
     if(p > 1.d-5)mu=sqrt(p) ! mu=cos theta - distance from disk center
	 mu_s(ix,iy)=mu
    enddo
  enddo

end

