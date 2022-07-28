subroutine model_spectrum(Jacobian,mode) ! calculate model spectrum in observed wavelengths by multiply LSD profile on line mask
use observed
use synthetic
use LSD_data
real(8) vp,prof,p,dv
real(8), parameter :: c=2.997925d5 
character(4), optional :: mode
real(8), optional :: Jacobian(:,:)

 model_prim=0.d0; model_sec=0.d0 ! Initialize model spectra arrays

 if(present(Jacobian))then
	if(mode == 'func')Jacobian=0.d0 ! Initialize Jacobian if need
 endif
 dv=vLSD_prim(2)-vLSD_prim(1)  ! Calculate model spectrum for primary
 do i=1,nlines_prim
  do k=1,nobs
   vp=(wave_obs(k)-w_mask_prim(i))/w_mask_prim(i)*c !  Corresponded for given wavelength velocity
   if(vp < vLSD_prim(1))cycle
   if(vp > vLSD_prim(nprim))cycle
   id=map1(vLSD_prim,rLSD_prim,nprim,vp,prof,1)
   model_prim(k)=model_prim(k)+prof*r_mask_prim(i)
   if(present(Jacobian))then
		if(mode == 'grad')then     ! Calculate Jacobian as function was linear interpolated
			ki=0
			do l=1,nprim
				if(vp <= vLSD_prim(l))then
					ki=l; exit
				endif
			enddo
			if(ki == 1)ki=2; if(ki == 0)ki=nprim
			p=(vp-vLSD_prim(ki-1))/dv
			Jacobian(k,ki-1)=Jacobian(k,ki-1)+r_mask_prim(i)*(1.d0-p)
			Jacobian(k,ki)=Jacobian(k,ki)+r_mask_prim(i)*p
		endif
   endif
  enddo
 enddo
 call correction(1) ! Correct model spectrum from synthetic spectrum

 dv=vLSD_sec(2)-vLSD_sec(1)  ! Calculate model spectrum for secondary
 do i=1,nlines_sec
  do k=1,nobs
   vp=(wave_obs(k)-w_mask_sec(i))/w_mask_sec(i)*c !  Corresponded for given wavelength velocity
   if(vp < vLSD_sec(1))cycle
   if(vp > vLSD_sec(nsec))cycle
   id=map1(vLSD_sec,rLSD_sec,nsec,vp,prof,1)
   model_sec(k)=model_sec(k)+prof*r_mask_sec(i)
   if(present(Jacobian))then
		if(mode == 'grad')then     ! Calculate Jacobian as function was linear interpolated
			ki=0
			do l=1,nsec
				if(vp <= vLSD_sec(l))then
					ki=l; exit
				endif
			enddo
			if(ki == 1)ki=2; if(ki == 0)ki=nsec
			p=(vp-vLSD_sec(l))/dv
			Jacobian(k,ki-1+nprim)=Jacobian(k,ki-1+nprim)+r_mask_sec(i)*(1.d0-p)
			Jacobian(k,ki+nprim)=Jacobian(k,ki+nprim)+r_mask_sec(i)*p
		endif
	endif
  enddo
 enddo

 call correction(2)
   	     
end

