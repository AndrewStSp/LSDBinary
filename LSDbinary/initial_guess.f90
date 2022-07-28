subroutine initial_guess(n) ! Choose the initial guess from synthetic spectra data
use LSD_data
use observed
use synthetic
use velocity
real(8) :: phi
 write(*,*)
 k=index(observed_names(n),slash,back=.true.) ! Link by file names in table
 write(*,*)trim(observed_names(n)(k+1:80))
 do i=1,number_spectra
  if(f_names(i) == observed_names(n)(k+1:80))then
   ni=i; exit
  endif 
 enddo
 BJD=JD(ni); phase=ph(ni)
 call radial_velocities(ni)

 phi = -100.d0
 if(nphases /= 0) then
	do i = 1, nphases
		if(phase == phases(i)) then
			phi = phase
			exit	
		endif
	enddo
	if(phi == -100.d0) write(*,*)'There is no such phase ...'
 endif

 if(phi == -100)write(*,'(3f10.3,a)')phase,Rv1,Rv2,'  (not eclipse)'
 if(phi /= -100)write(*,'(3f10.3,a)')phase,Rv1,Rv2,'  (eclipse)' 

 call read_corrandflux_file(phi)

 vLSD_prim=vLSD_syn1+Rv1; rLSD_prim=rLSD_syn1
 vLSD_sec=vLSD_syn2+Rv2; rLSD_sec=rLSD_syn2
 cp=cp_syn

end

subroutine radial_velocities(ni)
use LSD_data
use velocity
implicit real(8) (a-h,o-z)
real(8) M
data pi2/6.283185d0/

 select case(if_calc)

  case (0) ! Read from file
   Rv1=v1(ni); Rv2=v2(ni)
  case (1) ! Circlular orbit
   phase=ph(ni)
   t=-sin(pi2*phase)
   Rv1=K1*t+Vgamma
   Rv2=-K2*t+Vgamma

  case (2) ! Elliptical orbit
   phase=ph(ni)
   cos_omega=cos(omega*pi2/360.d0); sin_omega=sin(omega*pi2/360.d0); ecos_omega=e*cos_omega
   M=pi2*phase  ! Mean anomaly
   call solve_Kepler_equition(M,e,Ebig) ! Recieve eccentric anomaly E
   cos_E=cos(Ebig)
   cos_v=(cos_E-e)/(1.d0-e*cos_E) ! true anomaly
   sin_v=(sqrt(1.d0-e*e)*sin(Ebig))/(1.d0-e*cos_E)  
   Rv1=K1*(cos_v*cos_omega-sin_v*sin_omega+ecos_omega)+Vgamma
   omega2=omega-180.d0
   cos_omega=cos(omega2*pi2/360.d0); sin_omega=sin(omega2*pi2/360.d0); ecos_omega=e*cos_omega
   call solve_Kepler_equition(M,e,Ebig) ! Recieve eccentric anomaly E
   cos_E=cos(Ebig)
   cos_v=(cos_E-e)/(1.d0-e*cos_E) ! true anomaly
   sin_v=(sqrt(1.d0-e*e)*sin(Ebig))/(1.d0-e*cos_E)  
   Rv2=K2*(cos_v*cos_omega-sin_v*sin_omega+ecos_omega)+Vgamma
  
  end select
  
end

subroutine solve_Kepler_equition(M,e,Ebig)
implicit real(8) (a-h,o-z)
real(8) M
data maxiter/1000/
save

 if (e < 0.8d0)then  ! Initial guess
  Ebig=M
 else
  Ebig=M+sin(M)
 endif
 F=Ebig-e*sin(Ebig)-M
 iter=0
 do 
  Ebig=Ebig-F/(1.d0-(e*cos(Ebig)))
  F=Ebig-e*sin(Ebig)-M
  iter=iter+1
  if(abs(F) < 1.d-5.or.iter > maxiter)exit
 enddo
 if(iter > maxiter)write(*,*)'wrong Kepler equitiun solution!'

end
 
  
subroutine read_corrandflux_file(phi)
use synthetic
use LSD_data
implicit none
real(8) :: phi
integer :: ios, n, kkk, i
character(7) :: str

 1 continue
 if (phi == -100.d0) then
	open(2,file = trim(file_prim_corr),form = 'binary') !Load correction to model spectrum of primary
        read(2)nLSD_syn1,nsynth
	rewind (2)
	if(allocated(wcorr)) deallocate(wcorr,corr_prim,corr_sec,vLSD_syn1,rLSD_syn1)
	allocate (wcorr(nsynth), corr_prim(nsynth), corr_sec(nsynth),vLSD_syn1(nLSD_syn1),rLSD_syn1(nLSD_syn1),stat=ios)
	if(ios /= 0)stop 'Cannot allocate memory for correction'
	read(2)n,n,vLSD_syn1,rLSD_syn1,wcorr,corr_prim,cp_syn
	close (2)
 endif
 if (phi /= -100.d0) then
	kkk=index(file_prim_corr,'.',back=.true.)
    write(str,"(f7.4)")phi
	open(2,file = file_prim_corr(1:kkk-1)//'_'//trim(adjustl(str))//'.corr',status='old',form = 'binary',iostat=ios) !Load correction to model spectrum of primary
	if(ios/=0) then
		write(*,*)'file ',file_prim_corr(1:kkk-1)//'_'//trim(adjustl(str))//'.corr',' not found, will be calculated without taking into account the eclipse'
		phi = -100.d0
		goto 1
	endif
	read(2)nLSD_syn1,nsynth
	rewind (2)
	if(allocated(wcorr)) deallocate(wcorr,corr_prim,corr_sec,vLSD_syn1,rLSD_syn1)
	allocate (wcorr(nsynth), corr_prim(nsynth), corr_sec(nsynth),vLSD_syn1(nLSD_syn1),rLSD_syn1(nLSD_syn1),stat=ios)
	if(ios /= 0)stop 'Cannot allocate memory for correction'
	read(2)n,n,vLSD_syn1,rLSD_syn1,wcorr,corr_prim,cp_syn
	close (2)
 endif

 open(2,file = trim(file_sec_corr),form = 'binary') !Load correction to model spectrum of secondary
 read(2)nLSD_syn2
 rewind (2)
 if(allocated(vLSD_syn2))deallocate(vLSD_syn2,rLSD_syn2)
 allocate (vLSD_syn2(nLSD_syn2),rLSD_syn2(nLSD_syn2),stat=ios)
 read(2)n,n,vLSD_syn2,rLSD_syn2,wcorr,corr_sec
 close (2)
 nprim = nLSD_syn1; nsec = nLSD_syn2
 if(allocated(vLSD_prim))deallocate(vLSD_prim,rLSD_prim,vLSD_sec,rLSD_sec)
 allocate (vLSD_prim(nprim), rLSD_prim(nprim), vLSD_sec(nsec), rLSD_sec(nsec))

end subroutine read_corrandflux_file

