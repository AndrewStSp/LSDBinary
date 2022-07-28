subroutine load_observed 
use observed
real(8) w, r, w1, w2
integer(4) k, lll

 w1 = wfirst(1); w2 = wlast(number_regions)

 open(10,file=trim(adjustl(observed_spectrum)),status='old',iostat=ios)
 if(ios /= 0) then
  write(999,"('Observed spectrum ',a,' does not exist')") trim(adjustl(observed_spectrum)); stop
 endif
 lll=index(trim(adjustl(observed_spectrum)),'.rgs',back=.true.)
 nobs = 0  
 do          
  read(10,*,iostat=ios) w
  if(ios /= 0 .or. w > w2) exit
  if(w < w1) cycle
  do i = 1, number_regions
   if(w < wfirst(i) .or. w > wlast(i)) cycle
   nobs = nobs + 1
  enddo
 enddo
 allocate(wave_obs(nobs),R_obs(nobs)); rewind(10)

 nobs = 0
 do                  
  if(lll == 0)read(10,*,iostat = ios) w, r
  if(lll /= 0)read(10,*,iostat = ios) w, r, d, d, d, d
  if(ios /= 0 .or. w > w2) exit
  if(w < w1) cycle
  do i = 1, number_regions
   if(w < wfirst(i) .or. w > wlast(i)) cycle
   nobs = nobs + 1
   wave_obs(nobs) = w; R_obs(nobs) = 1.d0 - r
  enddo
 enddo
 close(10); wavelength_step = wave_obs(2) - wave_obs(1)

end



subroutine load_mask ! Load spectral lines masks (laboratory wavelength vs. line strengths)
use observed
use mask
character(4) elem_ID
real(8) w, r, w1, w2

 open(2,file=trim(adjustl(line_list)),status='old',iostat=ios)
 if(ios /= 0)then
  write(999,"('Line mask ',a,' does not exist')") trim(adjustl(line_list)); stop
 endif

 w1 = wfirst(1) - 2.d0; w2 = wlast(number_regions) + 2.d0; nlines = 0
 
 do 
  read(2,*,iostat=ios) w
  if(ios /= 0 .or. w > w2) exit
  if(w < w1) cycle
  do i = 1, number_regions
   if(w < (wfirst(i)-2.d0) .or. w > wlast(i)) cycle
   nlines = nlines + 1
  enddo
 enddo
 
 rewind(2); allocate(w_mask(nlines),r_mask(nlines),el_mask(nlines))
  
 m = 0
 do
  read(2,'(f10.4,1x,a4,f8.4)',iostat=ios) w, elem_ID, r
  if(ios /= 0 .or. w > w2) exit
  if(w < w1) cycle
  
  do i = 1, number_regions
   if(w < (wfirst(i)-2.d0) .or. w > wlast(i)) cycle
   m = m + 1
   w_mask(m) = w; r_mask(m) = 1.d0 - r; el_mask(m) = elem_ID
  enddo
 enddo
 close(2)

end
       