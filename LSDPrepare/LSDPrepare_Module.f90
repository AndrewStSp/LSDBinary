module init_prepare
    character(260)          atomic, molecular, model_prim, model_sec, conf_file
    character(80)           Star_name
    integer(4)              Resolve_power, nreg
    real(8), allocatable :: w1(:), w2(:)
    real                    Vsini_prim, Vsini_sec, Vturb_prim, Vturb_sec
    integer :: npoint_lsd
    real(8) :: min_prof
    contains
        character(260) function namef(textstring)
            character(*), intent(in) :: textstring
            integer                  :: k
            k=index(textstring,'\',back=.true.)
            if(k/=0) namef=' '//textstring(k+1:len(textstring))
        end function namef
end module init_prepare

module res
real,allocatable :: x(:), y(:)
integer :: chimg, check_lsd = 1, check_dat = 1, check_rgs = 1
end module res

module eclips
    real(8), allocatable :: phases(:)
    integer              :: nphases, if_eclips, covering
    real                 :: angle, ar, r2r1
end module eclips
