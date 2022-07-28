
program LSDPrerare

 use init_prepare
 use res
 use eclips

 npoint_lsd = 100; min_prof = 1.d-13
 open(999,file='log.txt')
 if_eclips = 0
 call read_init_conf
 call count_prepare

 close(999)
            
 stop

end program LSDPrerare

