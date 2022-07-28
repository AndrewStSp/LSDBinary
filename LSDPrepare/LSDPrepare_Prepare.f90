    subroutine count_prepare
      use init_prepare
      use LSD_data
      use res
!      use dflib
      use eclips
      character(200) string
      character(7) fn
      real(8) :: ff = -100

	  open(6,carriagecontrol = 'fortran')
      do i=1,2
! --------------------------------------------- SynthV --------------------------------------------------------------------------
        open(22,file='SynthV.config')
            write(22,'(a)')trim(atomic)
            write(22,'(a)')trim(molecular)
            if(i == 1)write(22,'((a))')trim(model_prim)
            if(i == 2)write(22,'((a))')trim(model_sec)
            write(22,*)w1(1),w2(nreg),' 0.01'
            if(i == 1)write(22,*)Vturb_prim
            if(i == 2)write(22,*)Vturb_sec
            if(i == 1)write(22,'((a))')trim(Star_name)//'_primary.syn'
            if(i == 2)write(22,'((a))')trim(Star_name)//'_secondary.syn'
            do m=1,4 
                write(22,*)'skip'
            enddo
        close(22)

        if(i == 1) write(999,*)'SynthV.config: ',trim(adjustl(atomic)),trim(adjustl(molecular)),&
                        w1(1),w2(nreg),' 0.01',trim(adjustl(model_prim)),Vturb_prim,trim(adjustl(Star_name))//'_primary.syn'

        if(i == 2) write(999,*)'SynthV.config: ',trim(adjustl(atomic)),trim(adjustl(molecular)),&
                        w1(1),w2(nreg),' 0.01',trim(adjustl(model_sec)),Vturb_sec,trim(adjustl(Star_name))//'_secondary.syn'

		write(6,*); string = ' '; write(string,"('./SynthV30 SynthV.config')")
!		write(6,*);	string = ' '; write(string,"('SynthV_main')")
		ios = system(trim(adjustl(string)))

        if(i == 1)then
            write(999,'((a))')'SynthV (primary star) ... OK'
			open(22,file=trim(adjustl(Star_name))//'_primary.syn'); close(22,status='delete')
        endif
        if(i == 2) then
            write(999,'((a))')'SynthV (secondary star) ... OK'
			open(22,file=trim(adjustl(Star_name))//'_secondary.syn'); close(22,status='delete')
        endif
		open(22,file='SynthV.config'); close(22,status='delete')

! --------------------------------------------- Convolve --------------------------------------------------------------------------
        open(22,file='convolve.conf')
            if(i == 1)write(22,'((a))')trim(Star_name)//'_primary.bsn'
            if(i == 2)write(22,'((a))')trim(Star_name)//'_secondary.bsn'
            if(i == 1)write(22,*)Vsini_prim
            if(i == 2)write(22,*)Vsini_sec
            write(22,*)Resolve_power
            if(i == 1)write(22,'((a))')trim(Star_name)//'_primary.rgs'
            if(i == 2)write(22,'((a))')trim(Star_name)//'_secondary.rgs'   
        close(22)

        if(i == 1) then
            write(999,*)'Convolve.conf: ',trim(adjustl(Star_name))//'_primary.bsn',Vsini_prim,Resolve_power,&
                        trim(adjustl(Star_name))//'_primary.rgs'
			write(6,'(1h+" Convolve (prim star) ... ",a4)')'wait'
        endif
        if(i == 2) then
            write(999,*)'Convolve.conf: ',trim(adjustl(Star_name))//'_secondary.bsn',Vsini_sec,Resolve_power,&
                        trim(adjustl(Star_name))//'_secondary.rgs'
			write(6,'(1h+" Convolve (sec star) ... ",a4)')'wait'
        endif

        call convolve

        if(i == 1) then
            write(999,'((a))')'Convolve (primary star) ... OK'
			write(6,'(1h+" Convolve (prim star) ... ",a4)')'OK  '; write(6,*)
        endif
        if(i == 2) then
            write(999,'((a))')'Convolve (secondary star) ... OK'
			write(6,'(1h+" Convolve (sec star) ... ",a4)')'OK  '; write(6,*)
			
        endif
! --------------------------------------------- LSDSynth --------------------------------------------------------------------------
        open(22,file='LSD.conf')
            if(i == 1)write(22,'((a))')trim(Star_name)//'_primary.rgs'
            if(i == 2)write(22,'((a))')trim(Star_name)//'_secondary.rgs'
            write(22,*)nreg
            do j=1,nreg
	           write(22,*)w1(j),w2(j)
            enddo
            if(i == 1)write(22,'((a))')trim(Star_name)//'_primary.lin'
            if(i == 2)write(22,'((a))')trim(Star_name)//'_secondary.lin'
        close(22)

        if(i == 1) then
            write(999,*)'LSD.conf: ',trim(adjustl(Star_name))//'_primary.rgs',nreg,'/',w1,'/','/',w2,'/',&
                        trim(Star_name)//'_primary.lin'
			write(6,'(1h+" LSDSynth (prim star) ... ",a4)')'wait'

            call LSDSynth(Vsini_prim)

            write(999,'((a))')'LSDSynth (prim star) ... OK'
			write(6,'(1h+" LSDSynth (prim star) ... ",a4)')'OK  '; write(*,*)
        endif
        if(i == 2) then
            write(999,*)'LSD.conf: ',trim(adjustl(Star_name))//'_secondary.rgs',nreg,'/',w1,'/','/',w2,'/',&
                        trim(Star_name)//'_secondary.lin'
			write(6,'(1h+" LSDSynth (sec star) ... ",a4)')'wait'

            call LSDSynth(Vsini_sec)

            write(999,'((a))')'LSDSynth (sec star) ... OK'
			write(6,'(1h+" LSDSynth (sec star) ... ",a4)')'OK  '; write(*,*)
			open(22,file=trim(adjustl(Star_name))//'_secondary.bsn'); close(22,status='delete')
        endif
! --------------------------------------------- end LSDSynth --------------------------------------------------------------------------

        if(if_eclips == 0)then        
            if(i == 1) then
				open(22,file=trim(adjustl(Star_name))//'_primary.bsn'); close(22,status='delete')
			endif
        endif

        deallocate(vLSD,rLSD)
		open(22,file=trim(adjustl(Star_name))//'_primary.lsd'); close(22,status='delete')
		open(22,file=trim(adjustl(Star_name))//'_secondary.lsd'); close(22,status='delete')
      enddo


      call flux_ratio(ff)
      if(if_eclips == 0)then
        write(*,*)'Complete'
        return
      endif
      do i=1,nphases
!		write(*,*)
        open(22,file='RME.config')
            write(22,*)r2r1 !'(f6.4)'
            write(22,*)Vsini_prim !'(f6.2)'
            write(22,*)angle
            write(22,*)ar
            write(22,*)Resolve_power !'(i6)'
            write(22,'((a))')trim(Star_name)//'_primary.bsn'
            write(22,'((a))')trim(Star_name)//'_secondary.rgs'
            write(22,*)phases(i)
        close(22)
        write(999,*)'RME.config: ',trim(Star_name)//'_primary.bsn',trim(Star_name)//'_secondary.rgs',&
                                                        r2r1,Vsini_prim,angle,ar,Resolve_power,phases(i)
		write(6,'(1h+" RME for phase:",f7.4," ... ",a4)')phases(i),'wait'

        call RME

        if(covering == 1) then  
			open(22,file=trim(Star_name)//'_primary_'//trim(adjustl(fn))//'.dat'); close(22,status='delete')
            string=''; write(string,'(a,f7.4)')' Eclipse star completely covers the eclipsed star in phase: ',phases(i)
			write(6,'(1h+"",a)')trim(string); write(*,*)
            write(999,'((a))')trim(string)
            cycle
        endif

        string=''; write(string,"('RME for phase: ',f7.4,' ... OK')")phases(i)
		write(999,'((a))')trim(string)
		write(6,'(1h+" RME for phase:",f7.4," ... ",a4)')phases(i),'OK  '; write(6,*)

        open(22,file='LSD.conf')
            write(fn,'(f7.4)')phases(i)
            write(22,'((a))')trim(Star_name)//'_primary_'//trim(adjustl(fn))//'.dat'
            write(22,*)nreg
            do j=1,nreg
	           write(22,*)w1(j),w2(j)
            enddo
            write(22,'((a))')trim(Star_name)//'_primary.lin'
        close(22)
        write(999,*)'LSD.conf (RME): ',trim(Star_name)//'_primary_'//trim(adjustl(fn))//'.dat',nreg,'/',w1,'/ /',w2,'/',&
                            trim(Star_name)//'_primary.lin'
		write(6,'(1h+" LSDSynth for eclips on phase:",f7.4," ... ",a4)')phases(i),'wait'


        call LSDSynth(Vsini_prim)

        string=''; write(string,"('LSDSynth for eclips on phase: ',f7.4,' ... OK')")phases(i)
        write(999,'((a))')trim(string)
		write(6,'(1h+" LSDSynth for eclips on phase:",f7.4," ... ",a4)')phases(i),'OK  '; write(*,*)


        call flux_ratio(phases(i))

        deallocate(vLSD,rLSD)
		open(22,file=trim(Star_name)//'_primary_'//trim(adjustl(fn))//'.dat'); close(22,status='delete')
		open(22,file=trim(Star_name)//'_primary_'//trim(adjustl(fn))//'.lsd'); close(22,status='delete')
      enddo

	  open(22,file=trim(adjustl(Star_name))//'_primary.bsn'); close(22,status='delete')
	  open(22,file=trim(Star_name)//'_secondary.rgs'); close(22,status='delete')
	  open(22,file=trim(Star_name)//'_primary.rgs'); close(22,status='delete')

      write(*,*)'Complete'

    end subroutine count_prepare


    subroutine read_init_conf
      use init_prepare
      use eclips
      character(2) :: s
	  integer ios
	  real aa
      open(1,file='LSDPrepare.conf',status='old',iostat=ios)
	  if(ios/=0)stop 'File LSDPrepare.conf not found'
      nreg=0; i=0; if_eclips=0; aa = 0.d0
      do
          read(1,'(a2)',iostat=ios)s
          if(ios /= 0)exit
          if(s(1:1) == '#')cycle
          backspace (1)
          select case (i)
              case (0)
                read(1,'((a))')Star_name; i=1
              case (1)
   				read(1,'((a))')atomic
				open(22,file=trim(atomic),status='old',iostat=ios)
				if(ios/=0) then
					write(*,*)'File ',trim(atomic),' not found'
					stop
				endif
				close(22); i=2
              case (2)
   	          read(1,'((a))')molecular; i=3 
              case (3)
   				read(1,*)Resolve_power
				if(Resolve_power < 0 .or. Resolve_power > 200000) then
					write(*,*)'Invalid Resolve power value (0 ... 200000)'
					stop
				endif
				i=4
              case (4)
                read(1,*)nreg
	           if(allocated(w1))deallocate(w1,w2); allocate(w1(nreg),w2(nreg))
                i=5
              case(5)
                do j=1,nreg
					read(1,*)w1(j),w2(j)
					aa = aa+(w2(j)-w1(j))
	            enddo
				if(aa > 1.d3) then
					write(999,*)'Interval selected too long. In total should not exceed 1000 A. Your length: ', aa
					write(*,*)'Interval selected too long. In total should not exceed 1000 A. Your length: ', aa
					stop
	            endif
				i=6
              case (6)
   				read(1,'((a))')model_prim
			  	open(22,file=trim(model_prim),status='old',iostat=ios)
				if(ios/=0) then
					write(*,*)'File ',trim(model_prim),' not found'
					stop
				endif
				close(22); i=7
              case (7)
   				read(1,*)Vsini_prim; i=8
              case (8)
   				read(1,*)Vturb_prim; i=9
              case (9)
   				read(1,'((a))')model_sec
			  	open(22,file=trim(model_sec),status='old',iostat=ios)
				if(ios/=0) then
					write(*,*)'File ',trim(model_sec),' not found'
					stop
				endif
				close(22); i=10
              case (10)
   				read(1,*)Vsini_sec; i=11
              case (11)
   				read(1,*)Vturb_sec; i=12
              case (12)    ! Eclipse
               read(1,*)s
			   if(s == 'E '.or.s == 'e ') then
			    if_eclips=1
				i=13
			   else
			    exit
			   endif
              case(13)
				read(1,*)r2r1; i=14
              case(14)
				read(1,*)angle; i=15
              case(15)
				read(1,*)ar; i=16
              case(16)
				read(1,*)nphases

	          if(allocated(phases))deallocate(phases); allocate(phases(nphases))
               i=17
              case(17)
				do k=1,nphases
	               read(1,*)phases(k)
                enddo
	            exit
          end select
      enddo
      close(1)
    end subroutine read_init_conf

     
  subroutine flux_ratio(phase) ! approximate continuum flux ratio of synthetic spectra
    use init_prepare
    real(8), allocatable :: w(:),f1(:),f2(:)
    real(8) cp(3),s(3),stat(10)
    real(8), intent(in) :: phase
    character(7) fn
    if_eclips=1
    if(phase < -10.d0)then
        if_eclips=0 ! not eclipse of primary
    endif
    if(if_eclips == 1) then
        write(fn,'(f7.4)')phase
        open(10,file=trim(Star_name)//'_primary_'//trim(adjustl(fn))//'.dat',iostat=ios)
    else
	   open(10,file=trim(Star_name)//'_primary.rgs',iostat=ios)
    endif
    n=0
    do
        read(10,*,iostat=ios)
        if(ios /= 0)exit
        n=n+1
    enddo
    rewind (10)
    allocate (w(n),f1(n),f2(n))

    do i=1,n
        read(10,*)w(i),d,d,f1(i)
    enddo
    close (10)
    open(10,file=trim(Star_name)//'_secondary.rgs',iostat=ios)
    do i=1,n
        read(10,*)d,d,d,f2(i)
    enddo
    close (10)
    f2=f2/f1
    w=log(w)
    call dRcurv(n,w,f2,2,cp,s,stat)
    if(if_eclips == 0) then
        open(10,file=trim(Star_name)//'_primary.corr', access='append', form = 'binary', iostat=ios)
	   if(ios /= 0) then
		write(999,*)'file = ',trim(Star_name)//'_primary.corr non exist'
		stop
	   endif
	   write(10)cp
    endif
    if(if_eclips /= 0) then
        write(fn,'(f7.4)')phase
        open(10,file=trim(Star_name)//'_primary_'//trim(adjustl(fn))//'.corr', access='append', form = 'binary', iostat=ios)
  	   if(ios /= 0) then
		write(999,*)'file = ',trim(Star_name)//'_primary_'//trim(adjustl(fn))//'.corr',' non exist'
		stop
	   endif
	   write(10)cp
    endif
    close (10)
    deallocate (w,f1,f2)
  end subroutine flux_ratio

    subroutine DRCURV(n,x,y,mp,d,s,st) ! Approximate y(x) by mp degree polinomial
      real(8) x(*),y(*),d(*),s(*),st(*)
      real(8),allocatable :: a(:),b(:),c(:),g(:),h(:),u(:),z(:),wg(:)
      integer(4),allocatable :: l(:)
      allocate (a(mp+1),b(mp+1),c(mp+1),g(mp+1),h(mp+1),u(2*n),l(mp+1),wg(n),z(n))
      lp=0
      wg=1.d0
      call Vc11ad(x,y,wg,z,n,a,b,c,g,h,l,mp,u,lp)
      call pe08ad(a,b,c,d,mp)
      deallocate (a,b,c,g,h,u,l,wg,z)
    end subroutine DRCURV

