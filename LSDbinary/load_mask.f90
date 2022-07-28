subroutine load_mask(n,file_name) ! Load spectral lines masks
use observed
use synthetic
character(4) elem_ID
character(80) file_name
real(8) w,r,w1,w2

 open(2,file=trim(file_name),status='old',iostat=ios)
 if(ios /= 0)then
  write(*,*)trim(file_name),'not exist'
  stop
 endif
 
 w1=wfirst(1)-2.d0; w2=wlast(number_regions)+2.d0
 m=0    ! check size of mask arrays and allocate it's
 do 
  read(2,*,iostat=ios)w
  if(w < w1)cycle
  if(w > w2)exit
  if(ios /= 0)exit
  do i=1,number_regions
   if(w < wfirst(i)-2.d0)cycle
   if(w > wlast(i))cycle
   m=m+1
  enddo
 enddo
 if(n == 1)then
  allocate(w_mask_prim(m),r_mask_prim(m),el_mask_prim(m))
  nlines_prim=m
 else
   allocate(w_mask_sec(m),r_mask_sec(m),el_mask_sec(m))
   nlines_sec=m
 endif
 rewind (2)

! Read mask data
 m=0
 do
  read(2,'(f10.4,1x,a4,f8.4)',iostat=ios)w,elem_ID,r
  if(w < w1)cycle
  if(w > w2)exit
  if(ios /= 0)exit
  do i=1,number_regions
   if(w < wfirst(i)-2.d0)cycle
   if(w > wlast(i))cycle
   m=m+1
   if(n == 1)then
    w_mask_prim(m)=w; r_mask_prim(m)=1.d0-r; el_mask_prim(m)=elem_ID(1:2)
   else
    w_mask_sec(m)=w; r_mask_sec(m)=1.d0-r; el_mask_sec(m)=elem_ID(1:2)
   endif
  enddo
 enddo

 close (2)

end
  