integer function map1 (xold, fold, nold, xnew, fnew, nnew) ! Interpolation function
 implicit real(8) (a-h,o-z)
 real(8) xold(*), fold(*), xnew(*), fnew(*)
 l  = 2
 ll = 0
 do k = 1, nnew
	do while ( xnew(k) >= xold(l) .and. l <= nold )
   	 l = l + 1
	enddo
	if ( l > nold ) then
		if ( l /= ll ) then
		
		   l  = amin0(nold,l)
		   c  = 0.d0
		   b  = (fold(l) - fold(l-1)) / (xold(l) - xold(l-1))
		   a  = fold(l) - xold(l) * b
		   ll = l
		
		endif
		
		fnew(k) = a + (b + c * xnew(k)) * xnew(k)
		cycle
	endif
	if ( l /= ll ) then
	    if ( l /= 2 .and. l /= 3 ) then
		   l1 = l - 1
		   if ( .not. ( l > ll+1 .or. ( l == 3 .or. l == 4 ) ) ) then
		
		      cbac = cfor
		      bbac = bfor
		      abac = afor
		    if ( l == nold ) then
		
			   c  = cbac
			   b  = bbac
			   a  = abac
			   ll = l
			   fnew(k) = a + (b + c * xnew(k)) * xnew(k)
			   cycle
		
		    endif
		else
			
		   l2 = l - 2
		   d  = (fold(l1) - fold(l2)) / (xold(l1) - xold(l2))
		   cbac = fold(l) / ((xold(l) - xold(l1)) * (xold(l) - xold(l2))) + &
				  (fold(l2) / (xold(l) - xold(l2)) - fold(l1) /             &
				  (xold(l) - xold(l1))) / (xold(l1) - xold(l2))
		   bbac = d - (xold(l1) + xold(l2)) * cbac
		   abac = fold(l2) - xold(l2) * d + xold(l1) * xold(l2) * cbac
		   if ( l >= nold ) then
		
			  c  = cbac
			  b  = bbac
			  a  = abac
			  ll = l
			  fnew(k) = a + (b + c * xnew(k)) * xnew(k)
			  cycle
		
		   endif
		endif
		d = (fold(l) - fold(l1)) / (xold(l) - xold(l1))
		cfor = fold(l+1) / ((xold(l+1) - xold(l)) * (xold(l+1) - xold(l1))) + &
			   (fold(l1) / (xold(l+1) - xold(l1)) - fold(l) /                 &
			   (xold(l+1) - xold(l))) / (xold(l) - xold(l1))
		bfor = d - (xold(l) + xold(l1)) * cfor
		afor = fold(l1) - xold(l1) * d + xold(l) * xold(l1) * cfor
		wt = 0.d0
		if ( abs(cfor) /= 0.d0 ) wt = abs(cfor) / ( abs(cfor) + abs(cbac) )
		
		   a  = afor + wt * (abac - afor)
		   b  = bfor + wt * (bbac - bfor)
		   c  = cfor + wt * (cbac - cfor)
		   ll = l
		   fnew(k) = a + (b + c * xnew(k)) * xnew(k)
		   cycle
		endif
		if ( l /= ll ) then
		
		   l  = amin0(nold,l)
		   c  = 0.d0
		   b  = (fold(l) - fold(l-1)) / (xold(l) - xold(l-1))
		   a  = fold(l) - xold(l) * b
		   ll = l
		
		endif
	endif
	fnew(k) = a + (b + c * xnew(k)) * xnew(k)
 enddo
 map1    = ll - 1
end function map1