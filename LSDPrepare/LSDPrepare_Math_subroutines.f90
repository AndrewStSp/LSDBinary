! Mathematical subroutines

pure subroutine MakeGaussMatrix(x1,x2,x,w,n)
implicit none
intent(in) :: x1,x2,n
intent(out) :: x,w
real(8), parameter :: pi=3.141525d0
real(8),parameter :: eps = 3.d-14 ! eps is the relative precision
integer(4) :: n,i,j,m
real(8)    :: x1,x2,x(n),w(n)
real(8)    :: p1,p2,p3,pp,xl,xm,z,z1

m = (n + 1)/2
xm = 0.5d0*(x2 + x1)
xl = 0.5d0*(x2 - x1)
do i=1,m
  z = cos(pi*(i-0.25d0)/(n+0.5d0))
  do while( .true. )
    p1 = 1.d0; p2 = 0.d0
    do j=1,n
      p3 = p2; p2 = p1; p1 = ((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
    enddo
    pp = n*(z*p1-p2)/(z*z-1.d0)
    z1 = z; z = z1-p1/pp            ! Newton's method.
    if( abs(z-z1)<=eps ) exit
  enddo
  x(i) = xm - xl*z                  ! Scale the root to the desired interval,
  x(n+1-i) = xm + xl*z              ! and put in its symmetric counterpart.
  w(i) = 2.d0*xl/((1.d0-z*z)*pp*pp) ! Compute the weight
  w(n+1-i) = w(i)                   ! and its symmetric counterpart.
enddo
end subroutine MakeGaussMatrix


subroutine integ ( x, f, fint, n, start)
  real(8) x(*), f(*), fint(*)
  real(8) A(n), B(n), C(n)
  real(8) start
  call parcoe ( f, x, A, B, C, n )
  fint(1) = start
  do i = 1, n - 1
     fint(i+1) = fint(i) + (A(i) + B(i) / 2.d0 * (x(i+1) + x(i)) + &
                 C(i) / 3.d0 * ((x(i+1) + x(i)) * x(i+1) + x(i) * x(i))) * (x(i+1) - x(i))
  enddo
end subroutine integ

subroutine parcoe ( f, x, A, B, C, n )
  real(8) f(*), x(*), A(*), B(*), C(*)
  real(8) wt, d
  C(1) = 0.d0
  B(1) = (f(2) - f(1)) / (x(2) - x(1))
  A(1) = f(1) - x(1) * B(1)
  C(n) = 0.d0
  B(n) = (f(n) - f(n-1)) / (x(n) - x(n-1))
  A(n) = f(n) - x(n) * B(n)
  if ( n > 2 ) then
	  do j = 2, n - 1
	     d    = (f(j) - f(j-1)) / (x(j) - x(j-1))
	     C(j) = f(j+1) / ((x(j+1) - x(j)) * (x(j+1) - x(j-1))) - f(j) / &
		        ((x(j) - x(j-1)) * (x(j+1) - x(j))) + f(j-1) /          &
				((x(j) - x(j-1)) * (x(j+1) - x(j-1)))
	     B(j) = d - (x(j) + x(j-1)) * C(j)
	     A(j) = f(j-1) - x(j-1) * d + x(j) * x(j-1) * C(j)
	  enddo
	  C(2) = 0.d0
	  B(2) = (f(3) - f(2)) / (x(3) - x(2))
	  A(2) = f(2) - x(2) * B(2)
	  C(3) = 0.d0
	  B(3) = (f(4) - f(3)) / (x(4) - x(3))
	  A(3) = f(3) - x(3) * B(3)
	  do j = 2, n - 1
	     if ( C(j) /= 0.d0 ) then
			wt   = abs(C(j+1)) / (abs(C(j+1)) + abs(C(j)))
			A(j) = A(j+1) + wt * (A(j) - A(j+1))
			B(j) = B(j+1) + wt * (B(j) - B(j+1))
			C(j) = C(j+1) + wt * (C(j) - C(j+1))
	     endif
	  enddo
	  A(n-1) = A(n)
	  B(n-1) = B(n)
	  C(n-1) = C(n)
  endif
end subroutine parcoe

subroutine linter ( xold, yold, nold, xnew, ynew, nnew )
  real(8) xold(*), yold(*), xnew(*), ynew(*)
  iold = 2
  do inew = 1, nnew
     do while( iold < nold .and. xold(iold) < xnew(inew) )
        iold = iold + 1
     enddo
     ynew(inew) = yold(iold-1) + (yold(iold) - yold(iold-1)) / &
	              (xold(iold) - xold(iold-1)) * (xnew(inew) - xold(iold-1))
  enddo
end subroutine linter

integer function map1 ( xold, fold, nold, xnew, fnew, nnew )
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




subroutine deriv ( x, f, dfdx, n )
  real(8) x(*), f(*), dfdx(*)
  real(8) s, scale, d1, d, Tan1, Tan
  dfdx(1) = (f(2) - F(1)) / (x(2) - x(1))
  dfdx(n) = (f(n) - f(n-1)) / (x(n) - x(n-1))
  if ( n > 2) then
	  s = abs( x(2) - x(1)) / (x(2) - x(1) )
	  do j = 2, n - 1
	     scale = dmax1(abs(f(j-1)),abs(f(j)),abs(f(j+1))) / abs(x(j))
	     if ( scale == 0.d0 ) scale = 1.d0
	     d1      = (f(j+1) - f(j)) / (x(j+1) - x(j)) / scale
	     d       = (f(j) - f(j-1)) / (x(j) - x(j-1)) / scale
	     Tan1    = d1 / (s * sqrt(1.d0 + d1**2) + 1.d0)
	     Tan     = d / (s * sqrt(1.d0 + d**2) + 1.d0)
	     dfdx(j) = (Tan1 + Tan) / (1.d0 - Tan1 * Tan) * scale
	  enddo
  endif
end subroutine deriv

real(8) function Voigt(v,a) ! Width5 subroutine, modified by T.Kipper
 implicit real(8) (a-h,o-z)
 dimension H0(41),H1(81),H2(41),Gaus20(20),Gaus10(10),wht10(10),&
           Gaus3(3),wht3(3),wht20(20)
 Data H0/ &
  1.0000000d0,0.9900500d0,0.9607890d0,0.9139310d0,0.8521440d0,0.7788010d0,&
  0.6976760d0,0.6126260d0,0.5272920d0,0.4448580d0,0.3678790d0,0.2981970d0,&
  0.2369280d0,0.1845200d0,0.1408580d0,0.1053990d0,0.0773050d0,0.0555760d0,&
  0.0391640d0,0.0270520d0,0.0183156d0,0.0121552d0,0.0079071d0,0.0050418d0,&
  0.0031511d0,0.0019305d0,0.0011592d0,0.0006823d0,0.0003937d0,0.0002226d0,&
  0.0001234d0,0.0000671d0,0.0000357d0,0.0000186d0,0.0000095d0,0.0000048d0,&
  0.0000024d0,0.0000011d0,0.0000005d0,0.0000002d0,0.0000001d0/
 Data H1/ &
  -1.1283800d0,-1.1059600d0,-1.0404800d0,-0.9370300d0,-0.8034600d0,&
  -0.6494500d0,-0.4855200d0,-0.3219200d0,-0.1677200d0,-0.0301200d0,&
   0.0859400d0, 0.1778900d0, 0.2453700d0, 0.2898100d0, 0.3139400d0,&
   0.3213000d0, 0.3157300d0, 0.3009400d0, 0.2802700d0, 0.2564800d0,&
   0.2317260d0, 0.2075280d0, 0.1848820d0, 0.1643410d0, 0.1461280d0,&
   0.1302360d0, 0.1165150d0, 0.1047390d0, 0.0946530d0, 0.0860050d0,&
   0.0785650d0, 0.0721290d0, 0.0665260d0, 0.0616150d0, 0.0572810d0,&
   0.0534300d0, 0.0499880d0, 0.0468940d0, 0.0440980d0, 0.0415610d0,&
   0.0392500d0, 0.0351950d0, 0.0317620d0, 0.0288240d0, 0.0262880d0,&
   0.0240810d0, 0.0221460d0, 0.0204410d0, 0.0189290d0, 0.0175820d0,&
   0.0163750d0, 0.0152910d0, 0.0143120d0, 0.0134260d0, 0.0126200d0,&
   0.0118860d0, 0.0112145d0, 0.0105990d0, 0.0100332d0, 0.0095119d0,&
   0.0090306d0, 0.0085852d0, 0.0081722d0, 0.0077885d0, 0.0074314d0,&
   0.0070985d0, 0.0067875d0, 0.0064967d0, 0.0062243d0, 0.0059688d0,&
   0.0057287d0, 0.0055030d0, 0.0052903d0, 0.0050898d0, 0.0049006d0,&
   0.0047217d0, 0.0045526d0, 0.0043924d0, 0.0042405d0, 0.0040964d0,&
   0.0039595d0/
  Data H2/ &
    1.0000000d0, 0.9702000d0, 0.8839000d0, 0.7494000d0, 0.5795000d0,&
    0.3894000d0, 0.1953000d0, 0.0123000d0,-0.1476000d0,-0.2758000d0,&
   -0.3679000d0,-0.4234000d0,-0.4454000d0,-0.4392000d0,-0.4113000d0,&
   -0.3689000d0,-0.3185000d0,-0.2657000d0,-0.2146000d0,-0.1683000d0,&
   -0.1282100d0,-0.0950500d0,-0.0686300d0,-0.0483000d0,-0.0331500d0,&
   -0.0222000d0,-0.0145100d0,-0.0092700d0,-0.0057800d0,-0.0035200d0,&
   -0.0021000d0,-0.0012200d0,-0.0007000d0,-0.0003900d0,-0.0002100d0,&
   -0.0001100d0,-0.0000600d0,-0.0000300d0,-0.0000100d0,-0.0000100d0,&
    0.0000000d0/
  Data Gaus20/ &
    0.05d0,0.15d0,0.25d0,0.35d0,0.45d0,0.6d0,0.75d0,0.9d0,1.05d0,1.2d0,&
    1.35d0,1.5d0,1.65d0,1.8d0,1.95d0,2.1d0,2.25d0,2.4d0,2.55d0,2.7d0/
  Data  wht20/ &
    0.0996d0,0.097634d0,0.0938184d0,0.0883723d0,0.100389d0,0.104549d0,&
    0.0855171d0,0.0668929d0,0.0500384d0,0.035796d0,0.0244891d0,0.0160216d0,&
    0.010024d0,0.00599769d0,0.00343183d0,0.00187789d0,0.000982679d0,&
    0.000491757d0,0.000235342d0,0.0000681124d0/
  Data Gaus10/ &
    0.245341d0,0.737474d0,1.234076d0,1.738538d0,2.254974d0,2.788806d0,&
    3.347855d0,3.944764d0,4.603682d0,5.38748d0/
  Data wht10/  &
    0.462244d0,0.2866755d0,0.1090172d0,0.02481052d0,0.003243773d0,&
    0.0002283386d0,0.7802556d-5,0.1086069d-6,0.4399341d-9,0.2229394d-12/
  Data Gaus3/ 0.436077d0,1.335849d0,2.35605d0/
  Data wht3/0.7246296d0,0.1570673d0,0.00453001d0/
 if (a <= 0.175)then
  v0=v*10.d0
  n=v0
   if(n >= 120)then
     vv=v*v
     Voigt=(0.56419+0.846/vv)/vv*a
    return
   endif
   if(n < 40)then
     v1=n
     n=n+1
     v2=v0-v1
     n1=n+1
     Voigt=v2*(H0(n1)-H0(n)+a*(H1(n1)-H1(n)+a*(H2(n1)-H2(n))))+&
           H0(n)+a*(H1(n)+a*H2(n))
    return
   endif
     n=n/2+20
     v1=(n-20)*2
     n=n+1
     v2=(v0-v1)/2.d0
     n1=n+1
     Voigt=A*((H1(n1)-H1(n))*v2+H1(n))
    return
 endif
 if(a <= 0.5d0.and.v <= 2.7d0)then ! Super 20 point quadrature
  Voigt=0.d0
   do i=1,20
    Voigt=Voigt+wht20(i)*A/3.1415925d0*(1.d0/(A*A+(v-Gaus20(i))**2)+&
                1.d0/(a*a+(v+Gaus20(i))**2))
   enddo
  return
 endif
 if((a <= 0.5d0.and.v <= 7.d0).or.(a <= 1.d0.and.v <= 2.7d0))then ! 10 point
  Voigt=0.d0
   do i=1,10
    Voigt=Voigt+wht10(i)*a/3.1415925d0*(1.d0/(a*a+(v-Gaus10(i))**2)+&
                1.d0/(A*A+(v+Gaus10(i))**2))
   enddo
  return
 endif
 if((a <= 0.5d0.and.v <=100.d0).or.(a <= 1.d0.and.v <= 40.d0).or.&
    (a <= 15.d0.and.v <= 20.d0))then ! 3 point Gaussian quadrature
  Voigt=0.d0
   do i=1,3
    Voigt=Voigt+wht3(i)*A/3.1415925d0*(1.d0/(a*a+(v-Gaus3(i))**2)+&
                1.d0/(A*A+(v+Gaus3(i))**2))
   enddo
  return
 endif
   Voigt=a/1.77245d0/(a*a+v*v)  ! Lorentzian
end


real(8) function expi (n, x)
 implicit none
  real(8), parameter :: A0 = -44178.5471728217d0
  real(8), parameter :: A1 =  57721.7247139444d0
  real(8), parameter :: A2 =  9938.31388962037d0
  real(8), parameter :: A3 =  1842.11088668000d0
  real(8), parameter :: A4 =  101.093806161906d0
  real(8), parameter :: A5 =  5.03416184097568d0
  real(8), parameter :: B0 =  76537.3323337614d0
  real(8), parameter :: B1 =  32597.1881290275d0
  real(8), parameter :: B2 =  6106.10794245759d0
  real(8), parameter :: B3 =  635.419418378382d0
  real(8), parameter :: B4 =  37.2298352833327d0
  real(8), parameter :: C0 =  4.65627107975096d-7
  real(8), parameter :: C1 =  0.999979577051595d0
  real(8), parameter :: C2 =  9.04161556946329d0
  real(8), parameter :: C3 =  24.3784088791317d0
  real(8), parameter :: C4 =  23.0192559391333d0
  real(8), parameter :: C5 =  6.90522522784444d0
  real(8), parameter :: C6 =  0.430967839469389d0
  real(8), parameter :: D1 =  10.0411643829054d0
  real(8), parameter :: D2 =  32.4264210695138d0
  real(8), parameter :: D3 =  41.2807841891424d0
  real(8), parameter :: D4 =  20.4494785013794d0
  real(8), parameter :: D5 =  3.31909213593302d0
  real(8), parameter :: D6 =  0.103400130404874d0
  real(8), parameter :: E0 = -0.999999999998447d0
  real(8), parameter :: E1 = -26.6271060431811d0
  real(8), parameter :: E2 = -241.055827097015d0
  real(8), parameter :: E3 = -895.927957772937d0
  real(8), parameter :: E4 = -1298.85688746484d0
  real(8), parameter :: E5 = -545.374158883133d0
  real(8), parameter :: E6 = -5.66575206533869d0
  real(8), parameter :: F1 = 28.6271060422192d0
  real(8), parameter :: F2 = 292.310039388533d0
  real(8), parameter :: F3 = 1332.78537748257d0
  real(8), parameter :: F4 = 2777.61949509163d0
  real(8), parameter :: F5 = 2404.01713225909d0
  real(8), parameter :: F6 = 631.657483280800d0
  real(8) :: x1 = -1.0d20
  real(8)    x, ex, ex1
  integer n, i


  if ( x /= x1 ) then

    ex = exp( - x )
    x1 = x
    if ( x > 4.d0 ) then
      ex1 = (ex + ex * (E0 + (E1 + (E2 + (E3 + (E4 + (E5 + E6 / x) / x) / x) /&
             x) / x) / x) / &
             (x + F1 + (F2 + (F3 + (F4 + (F5 + F6 / x) / x) / x) / x) / x)) / x
    else if ( x > 1.d0 ) then
      ex1 = ex * (C6 + (C5 + (C4 + (C3 + (C2 + (C1 + C0 * x) * x) * x) * x) * &
            x) * x) / (D6 + (D5 + (D4 + (D3 + (D2 + (D1 + x) * x) * x) * x) * &
            x) * x)
    else if ( x > 0.d0 ) then
      ex1 = (A0 + (A1 + (A2 + (A3 + (A4 + A5 * x) * x) * x) * x) * x) / &
            (B0 + (B1 + (B2 + (B3 + (B4 + x) * x) * x) * x) * x) - log(x)
    else
      ex1 = 0.d0
    endif

  endif

  expi = ex1

  if ( n > 1 ) then

    do i = 1, n - 1
      expi = (ex - x * expi) / dble( i )	
    enddo

  endif
end function expi



subroutine solvit(A,n,b,ipivot)
 implicit real*8 (a-h,o-z)
 dimension A(n,n),b(n),ipivot(n)
  n1=n-1
   do i=1,n1
    m=i
    i1=i+1
    do  k=i1,n
     if(abs(A(k,i)) > abs(A(m,i)))m=k
    enddo
   ipivot(i)=m
   if(m /= i)then
    do k=i1,n
     t=A(i,k)
     A(i,k)=A(m,k)
     A(m,k)=t
    enddo
   endif
   pivot=1.d0/A(m,i)
   A(m,i)=A(i,i)
   a(i,i)=pivot
    do k=i1,n
     a(k,i)=a(k,i)*pivot
    enddo
    do j=i1,n
     c=A(i,j)
     if(c /= 0.d0)then
      do  k=i1,n
       A(k,j)=A(k,j)-A(k,i)*c
      enddo
     endif
    enddo
  enddo
   a(n,n)=1.d0/a(n,n)
  do i=1,n1
   m=ipivot(i)
   if(m /= i)then
    t=b(m)
    b(m)=b(i)
    b(i)=t
   endif
   c=b(i)
   i1=i+1
    do k=i1,n
     b(k)=b(k)-A(k,i)*c
    enddo
  enddo
  j1=n
  do i=1,n1
   j=j1
   j1=j1-1
   b(j)=b(j)*A(j,j)
   c=b(j)
    do k=1,j1
     b(k)=b(k)-A(k,j)*c
    enddo
  enddo
   b(1)=b(1)*A(1,1)
end 




subroutine mini(xopt,obs,n,nwl,fun)

real(8) xopt(n), obs(nwl), weight(nwl)
real(8) ftoler, freduc, fvalue, fopt, dev, oldfun, d, regul, cormax, corlim, corscl, di, qlambda, qscale
logical failure, done
integer info, wrksiz, iterf, maxitr, iterc, nspec
parameter (ftoler = 0.0005d0, freduc = 0.9995d0, corlim = 3.d0)
real(8), allocatable :: x(:), prf(:), work(:), sigc(:,:), diag(:), rhsold(:), deriva(:,:), derivc(:,:), rhs(:)
integer, allocatable :: ipiv(:)

allocate (x(n), prf(nwl), ipiv(n), work(45*n), sigc(n,n), diag(n), rhsold(n), deriva(nwl,n), derivc(n,n), rhs(n), stat=ios)
if(ios /= 0)stop 'Cannot allocate arrays in mini subroutine'
work = 0.d0; iterc = 0; deriva = 0.d0; wrksiz = 45*n; maxitr = 50; ncalls = 0; failure = .false.; qlambda = 10.0
qscale  = 10.0; fopt = 1.d30; iterf = iterc; weight = 1.d0
do 
 iterc = iterc + 1
 call dcopy(n,xopt,1,x,1)
 oldfun = fopt
 call fun(x,prf,deriva,n,nwl,'grad')
 do i = 1, n
  do iwl = 1, nwl
   deriva(iwl,i) = deriva(iwl,i) * sqrt(weight(iwl) / nwl)
  enddo
 enddo 
 do i = 1, nwl
 enddo
 call deviat(obs,prf,weight,nwl,d,dev)
 fvalue = d
 if(fvalue < fopt) then
  call dcopy(n,x,1,xopt,1)
  fopt = fvalue
 endif
 call dsyrk('u','t',n,nwl,1.d0,deriva,nwl,0.d0,derivc,n)
 do i1 = 2, n
  do i2 = 1, i1-1
   derivc(i1,i2) = derivc(i2,i1)
  enddo
 enddo
 do i = 1, n
  di = 0.d0
  do iwl = 1, nwl
   di = di + sqrt(weight(iwl)) * (obs(iwl) - prf(iwl)) * deriva(iwl,i)
  enddo
  rhs(i) = di / sqrt(dble(nwl))
 enddo
 call dcopy(n,derivc,n+1,diag,1)
 call dcopy(n,rhs,1,rhsold,1)
 done = .false.
 ifail = 0
 itry = 0
 oldfun = fvalue
 do
  call dscal(n,1.d0+qlambda,derivc,n+1)
  ifail = ifail + 1
  itry = itry + 1 
  call dsysv('u',n,1,derivc,n,ipiv,rhs,n,work,wrksiz,info)
  if(info /= 0) failure = .true.
  if(failure) then
   write(999,'('' Degenerate matrix DERIVC.'','' Cannot solve SLE.'',T95,''*** (MINI  ) ***'')')
   write(999,*) 'INFO=', info  
   deallocate (x, prf, ipiv, work, sigc, diag, rhsold, deriva, derivc, rhs)
   stop
  endif
  cormax = 0.d0
  do i = 1, n
   if(cormax < abs(rhs(i))) cormax = abs(rhs(i))
  enddo
  corscl = 1.d0
  if(cormax > corlim) corscl = corlim / cormax
  do i = 1, n
   x(i) = x(i) + rhs(i) * corscl
  enddo
  call fun(x,prf,deriva,n,nwl,'func')
  call deviat(obs,prf,weight,nwl,d,dev)
  fvalue = d; dev = sqrt(dev)
  if(itry == 2 .and. fvalue > oldfun*freduc) then
   if(qscale > 1.) then
    qscale = 0.1
   else
    qscale = 10.
   endif
   qlambda = qlambda * qscale * qscale
   call dcopy(n,xopt,1,x,1)
   fvalue = fopt
   done = .false.
  elseif(itry > 2 .and. fvalue > oldfun*freduc) then
   qlambda = qlambda / qscale
   done = .true.   
  else
   qlambda = qlambda * qscale
   if(qlambda > 1.e10) then
    qlambda = 1.e9
	done = .true.
   elseif(qlambda < 1.e-10) then
    qlambda = 1.e-9
	done = .true.
   else
    done = .false.
   endif
  endif
  if(fvalue < fopt*(1.d0-ftoler)) then
   call deviat(obs,prf,weight,nwl,di,dev)
   dev = sqrt(dev)
   call dcopy(n,x,1,xopt,1)
   fopt = fvalue
   ifail = 0
  endif
  do i1 = 2, n
   do i2 = 1, i1-1
    derivc(i2,i1) = derivc(i1,i2)
   enddo
  enddo
  call dcopy(n,diag,1,derivc,n+1)
  call dcopy(n,rhsold,1,rhs,1)
  if(.not. done .or. itry <= 2 .or. fvalue < oldfun*freduc) then
   oldfun = fvalue
   cycle
  else
   exit ! From internal loop
  endif
 enddo
 write(999,"(' ====  ITER:',I3,' MEAN DEV.=',F7.3,'% FUNCT.=',G12.5,' FOPT=',G12.5)") iterc, dev*100, fvalue, fopt
 if(ifail < itry .or. (iterc - iterf) < 2) then
  call dcopy(n,xopt,1,x,1)
  fvalue = fopt
  oldfun = fopt 
  if(iterc >= maxitr) then
   write(999,'(''The limit of iteration number has been reached'',t95,''*** (mini  ) ***'')')
   call deviat(obs,prf,weight,nwl,di,dev)     
   dev = sqrt(dev)
   exit ! From main loop
  endif     
 else
  exit ! From main loop
 endif 
enddo ! End of main loop     
call dcopy(n,xopt,1,x,1)
deallocate (x, prf, ipiv, work, sigc, diag, rhsold, deriva, derivc, rhs)
end

!------------------------------------------------------------------------------


subroutine deviat(obs,prf,weight,nwl,dev1,dev2)
real(8) prf(nwl), obs(nwl), weight(nwl), dev1, dev2
dev1 = 0.d0; dev2 = 0.d0
do iwl = 1,nwl
 dev1 = dev1+weight(iwl)*(prf(iwl)-obs(iwl))**2
 dev2 = dev2+(prf(iwl)-obs(iwl))**2
enddo
dev1 = dev1/nwl; dev2 = dev2/nwl
end

subroutine pe08ad (a,b,c,d,m) 
 double precision a,b,c,cd,cm,d,w1,w2 
 dimension cd(200),a(*),b(*),c(*),d(*) 
      if(m)99,1,2 
    1 d(1)=c(1) 
      return 
    2 if(m-1)99,3,4 
    3 d(1)=c(1)-a(1)*c(2) 
      d(2)=c(2) 
      return 
    4 cm=c(m+1) 
      cd(1)=c(m)-a(m)*cm 
      d(1)=c(m-1)-b(m)*cm 
      if (m-2)99,90,5 
    5 m1=m-1 
      w1=cd(1) 
      cd(1)=d(1)-a(m-1)*w1 
      cd(2)=w1-cm*a(m-1) 
      d(1)=c(m-2)-w1*b(m-1) 
      d(2)=-cm*b(m-1) 
      if (m-3) 99,90,6 
    6 do 60 i=3,m1 
      ir=m-i+1 
      w1=cd(1) 
      cd(1)=d(1)-w1*a(ir) 
      cd(i)=cd(i-1)-cm*a(ir) 
      d(1)=c(ir-1)-w1*b(ir) 
      d(i)=-cm*b(ir) 
      i1=i-1 
      do 70 j=2,i1 
      w2=cd(j) 
      cd(j)=w1+d(j)-w2*a(ir) 
      d(j)=-w2*b(ir) 
   70 w1=w2 
   60 continue 
   90 d(1)=d(1)-a(1)*cd(1) 
      d(m)=cd(m-1)-a(1)*cm 
      d(m+1)=cm 
      if(m-2) 99,99,91 
   91 do 92 j=2,m1 
   92 d(j)=cd(j-1)+d(j)-a(1)*cd(j) 
   99 return 
end                                           
                                                                        
                                                                        
                                                                        
subroutine VC11AD(X,Y,W,Z,N,A,B,C,G,H,L,M,U,LP) 
 integer LP,M,N 
 double precision A(*),B(*),C(*),G(*),H(*),U(*),W(*),X(*),Y(*),Z(*) 
 integer L(*) 
 double precision EJ,FJ,GJ,H1,W1,W2,WS,YY 
 integer I,II,IV,J,J1,L2,M1 
 intrinsic DSIGN,DSQRT,FLOAT,IDINT 
      EJ = 0.0D0 
      FJ = 0.0D0 
      GJ = 0.0D0 
      do 1 I = 1,N 
        IV = I + N 
        U(I) = dsqrt(W(I)) 
        Z(I) = Y(I)*U(I) 
        EJ = EJ + X(I)*U(I)**2 
        FJ = FJ + U(I)*Z(I) 
        GJ = GJ + U(I)**2 
        U(IV) = 0.0D0 
    1 enddo 
      A(1) = EJ/GJ 
      B(1) = 0.0D0 
      C(1) = FJ/GJ 
      G(1) = 1.0D0/GJ 
      do 2 J = 1,M 
        EJ = 0.0D0 
        FJ = 0.0D0 
        GJ = 0.0D0 
        H1 = 0.0D0 
        L2 = 0 
        do 3 I = 1,N 
          Z(I) = Z(I) - C(J)*U(I) 
          IV = I + N 
          WS = (X(I)-A(J))*U(I) - B(J)*U(IV) 
          U(IV) = U(I) 
          U(I) = WS 
          EJ = EJ + X(I)*U(I)**2 
          FJ = FJ + U(I)*Z(I) 
          GJ = GJ + U(I)**2 
          H1 = H1 + Z(I)**2 
          if (1-I) 13,3,3 
   13     L2 = L2 + IDINT((DSIGN(1.0D0,-Z(I)*Z(I-1))+1.0D0)/2.0D0) 
    3   continue 
        A(J+1) = EJ/GJ 
        B(J+1) = GJ*G(J) 
        C(J+1) = FJ/GJ 
        G(J+1) = 1.0D0/GJ 
        H(J) = H1 
        L(J) = L2 
    2 enddo 
      M1 = M + 1 
      if (LP.GT.0) write (LP,FMT=20) 
   20 format (1H1,8X,4HX(I),12X,4HY(I),12X,4HZ(I),/,/) 
      H1 = 0.0D0 
      L2 = 0 
      do 10 I = 1,N 
        W2 = 0.0D0 
        YY = C(M+1) 
        do 101 J = 1,M 
          J1 = M - J + 1 
          W1 = W2 
          W2 = YY 
          YY = C(J1) + (X(I)-A(J1))*W2 - B(J1+1)*W1 
  101   continue 
        Z(I) = YY 
        U(I) = (Z(I)-Y(I)) 
        H1 = H1 + U(I)*U(I)*W(I) 
        if (1-I) 111,110,110 
  111   L2 = L2 + IDINT((DSIGN(1.0D0,-U(I)*U(I-1))+1.0D0)/2.0D0) 
  110   if (LP.GT.0) write (LP,FMT=11) X(I),Y(I),Z(I) 
   10 enddo 
      H(M1) = H1 
      L(M1) = L2 
      H1 = H1/FLOAT(N-M1) 
      do 21 J = 1,M1 
        G(J) = G(J)*H1 
   21 enddo 
      if (LP.LE.0) GO TO 40 
      write (LP,FMT=7) 
    7 format (1H1,61X,8HVARIANCE,10X,9HRESIDUALS) 
      write (LP,FMT=6) 
    6 format (7X,1HJ,7X,4HA(J),12X,4HB(J),12X,4HC(J),11X,7HOF C(J),6X,  &
     &       14HSUM OF SQUARES,3X,12HSIGN CHANGES)                      
      M1 = M + 1 
      DO 24 I = 1,M1 
        II = I - 1 
   24 write (LP,FMT=4) II,A(I),B(I),C(I),G(I),H(I),L(I) 
   11 format (3E16.6) 
   40 return 
    4 format (3X,I5,5E16.6,6X,I5) 
end                                           
                                                                        
                                                                        
SUBROUTINE VB05BD(M,N,XD,YD,WD,RD,XN,FN,GN,IPRINT,W) 
 INTEGER IPRINT,M,N 
 DOUBLE PRECISION FN(*),GN(*),RD(*),W(*),WD(*),XD(*),XN(*),YD(*) 
 DOUBLE PRECISION A,AA,AM,B,C,D 
 REAL HPPR 
 INTEGER I,II,IL,IS,IU,J,JJ,JU,K,KA,KB,KC,KK,KL,ML,MM,MU,NINC,NN,NU 
 DOUBLE PRECISION S(5) 
 INTRINSIC DABS,DSIGN,DSQRT,DBLE,IFIX,MAX0,SNGL 
      ML = N + N + 13 
      W(ML-6) = XN(1) + XN(1) - XN(2) 
      K = ML 
      DO 1 I = 1,M 
        IF (WD(I)) 2,1,2 
    2   IF (XD(I)-W(K-6)) 3,3,4 
    3   A = DSQRT(W(K-2)**2+WD(I)**2) 
        W(K-1) = (W(K-2)*W(K-1)+YD(I)*WD(I)**2)/A 
        W(K-2) = A 
        GO TO 1 
    4   W(K) = XD(I) 
        W(K+5) = YD(I)*WD(I) 
        W(K+4) = WD(I) 
        K = K + 6 
    1 END DO 
      MM = (K-ML)/6 
      MU = K - 1 
      IF (MM-4) 5,6,6 
    5 WRITE (6,FMT=7) MM 
    7 FORMAT (/,/,5X,21HVB05BD THERE ARE ONLY,I2,18H INDEPENDENT DATA,, &
     &      57H SO A POLYNOMIAL OF DEGREE LESS THAN THREE WILL BE FITTED&
     &       ,/,/)                                                      
      W(1) = XN(1) 
      W(2) = XN(N) 
      A = 0.0D0 
      IF (MM) 8,8,9 
    9 A = W(ML+5)/W(ML+4) 
      GO TO (8,10,11) MM 
    8 DO 12 I = 1,N 
        FN(I) = A 
        GN(I) = 0.0D0 
   12 END DO 
      GO TO 13 
   10 B = (W(ML+11)/W(ML+10)-A)/ (W(ML+6)-W(ML)) 
      DO 14 I = 1,N 
        FN(I) = A + B* (XN(I)-W(ML)) 
        GN(I) = B 
   14 END DO 
      GO TO 13 
   11 B = W(ML+11)/W(ML+10) 
      A = (B-A)/ (W(ML+6)-W(ML)) 
      D = (W(ML+17)/W(ML+16)-B)/ (W(ML+12)-W(ML+6)) 
      C = (A* (W(ML+12)-W(ML+6))+D* (W(ML+6)-W(ML)))/ (W(ML+12)-W(ML)) 
      D = (D-A)/ (W(ML+12)-W(ML)) 
      DO 15 I = 1,N 
        FN(I) = ((XN(I)-W(ML+6))*D+C)* (XN(I)-W(ML+6)) + B 
        GN(I) = C + 2.0D0*D* (XN(I)-W(ML+6)) 
   15 END DO 
      GO TO 13 
    6 W(MU+1) = XN(N) + XN(N) - XN(N-1) 
      NN = 4 
      J = ML 
      K = 3 
      NINC = N + 6 
      W(NN) = XN(1) 
      W(NINC+4) = DBLE(ML) - 5.5D0 
      DO 17 I = 2,N 
        IF (K-3) 18,18,19 
   19   WRITE (6,FMT=20) W(NN) 
   20   FORMAT (/,/,5X,                                                 &
     &         49HVB05BD THERE ARE TOO MANY KNOTS SO THE ONE AT X =,    &
     &         E13.4,17H HAS BEEN DELETED,/,/)                          
        GO TO 21 
   18   NN = NN + 1 
        K = K + 1 
   21   IF (XN(I)-W(J)) 22,22,23 
   23   K = K - 1 
        J = J + 6 
        GO TO 21 
   22   W(NN) = XN(I) 
        NU = NN + NINC 
        W(NU) = DBLE(J) - 5.5D0 
        K = MAX0(K,0) 
   17 END DO 
      K = K - (MU+1-J)/6 
   16 IF (K) 24,24,25 
   25 K = K - 1 
      NN = NN - 1 
      WRITE (6,FMT=20) W(NN) 
      GO TO 16 
   24 W(NN) = XN(N) 
      NU = NN + NINC 
      W(NU) = DBLE(MU) - 4.5D0 
      DO 26 I = 1,3 
        J = NN + I 
        W(J) = W(J-1) + W(J-1) - W(J-2) 
        NU = J + NINC 
        W(NU) = W(NU-1) 
        J = 4 - I 
        W(J) = W(J+1) + W(J+1) - W(J+2) 
        NU = J + NINC 
        W(NU) = W(NU+1) 
   26 END DO 
      IS = 1 
      I = 3 
      GO TO 28 
   29 NU = I + NINC 
      HPPR = SNGL(W(NU-2)) 
      K = 6 + IFIX(HPPR) 
   30 IF (DBLE(K)-W(NU-1)) 31,31,32 
   31 W(K+4) = W(K+4)*A* (W(K)-W(I-2))**3 
      IF (DABS(W(K+4)).LE.1.0D-37) W(K+4) = 0.0D0 
      K = K + 6 
      GO TO 30 
   32 IF (DBLE(K)-W(NU)) 33,33,34 
   33 W(K+3) = W(K+4)* (A* (W(K)-W(I-2))**3+B* (W(K)-W(I-1))**3) 
      IF (DABS(W(K+3)).LE.1.0D-37) W(K+3) = 0.0D0 
      K = K + 6 
      GO TO 32 
   34 IF (DBLE(K)-W(NU+1)) 35,35,36 
   35 W(K+2) = W(K+4)* (C* (W(I+1)-W(K))**3+D* (W(I+2)-W(K))**3) 
      IF (DABS(W(K+2)).LE.1.0D-37) W(K+2) = 0.0D0 
      K = K + 6 
      GO TO 34 
   36 IF (DBLE(K)-W(NU+2)) 37,37,38 
   37 W(K+1) = W(K+4)*D* (W(I+2)-W(K))**3 
      IF (DABS(W(K+1)).LE.1.0D-37) W(K+1) = 0.0D0 
      K = K + 6 
      GO TO 36 
   38 IF (I-NN) 39,39,27 
   39 I = I + 1 
      GO TO 28 
   27 II = ML - 6 
      J = NINC 
   40 J = J + 1 
      IF (NINC+NN-J) 66,66,65 
   65 II = II + 6 
      K = 4 
      IF (W(J)-DBLE(II)) 41,41,42 
   41 JJ = J - K 
      IF (W(JJ+5)-DBLE(II)) 43,43,44 
   43 K = K - 1 
      GO TO 41 
   42 W(II+4) = 0.0D0 
      JJ = J - K 
   44 HPPR = SNGL(W(JJ+4)) 
      IL = 6 + MAX0(II,IFIX(HPPR)) 
      HPPR = SNGL(W(JJ+5)) 
      IU = IFIX(HPPR) 
      IF (IU-IL) 67,68,68 
   67 K = K - 1 
      IF (K) 40,40,69 
   69 W(II+1) = W(II+2) 
      W(II+2) = W(II+3) 
      W(II+3) = W(II+4) 
      GO TO 42 
   68 KB = II + K 
      DO 70 KK = K,5 
        KA = II + KK 
        S(KK) = W(KA)*W(KB) 
   70 END DO 
      DO 71 I = IL,IU,6 
        KC = I + K 
        DO 72 KK = K,5 
          KA = I + KK 
          S(KK) = S(KK) + W(KA)*W(KC) 
   72   CONTINUE 
   71 END DO 
      AA = DSIGN(DSQRT(S(K)),W(KB)) 
      AM = S(K) + AA*W(KB) 
      KL = K + 1 
      DO 73 KK = KL,5 
        KA = II + KK 
        S(KK) = (S(KK)+AA*W(KA))/AM 
        W(KA) = W(KA) - S(KK)* (W(KB)+AA) 
   73 END DO 
      W(KB) = -AA 
      DO 74 I = IL,IU,6 
        KC = I + K 
        DO 75 KK = KL,5 
          KA = I + KK 
          W(KA) = W(KA) - S(KK)*W(KC) 
   75   CONTINUE 
   74 END DO 
      GO TO 67 
   66 JU = J - 2 
   45 W(JU+1) = (W(II+5)-W(JU+2)*W(II+2)-W(JU+3)*W(II+3)-               &
     &          W(JU+4)*W(II+4))/W(II+1)                                
      II = II - 6 
      JU = JU - 1 
      IF (JU-NINC) 46,45,45 
   46 IS = 2 
      I = 3 
      K = ML 
      J = NINC 
      GO TO 28 
   47 J = J + 1 
      W(K) = W(K) + W(J)*A* (W(I-1)-W(I-2))**3 
      W(K+1) = 3.D0* (W(K+1)+W(J)*A* (W(I-1)-W(I-2))**2) 
      W(K+2) = W(K+2) + W(J)* (A* (W(I)-W(I-2))**3+B* (W(I)-W(I-1))**3) 
      W(K+3) = W(K+3) + W(J)* (A* (W(I)-W(I-2))**2+B* (W(I)-W(I-1))**2) 
      W(K+4) = W(J)*D* (W(I+2)-W(I+1))**3 
      W(K+5) = -W(J)*D* (W(I+2)-W(I+1))**2 
      W(I-2) = W(I+1) 
      IF (NN-I) 48,49,49 
   49 I = I + 1 
      K = K + 2 
      GO TO 28 
   48 J = 1 
      K = ML + 4 
      DO 50 I = 1,N 
        IF (W(J)-XN(I)) 51,51,52 
   51   FN(I) = W(K) 
        GN(I) = W(K+1) 
        J = J + 1 
        K = K + 2 
        GO TO 50 
   52   D = (XN(I)-W(J-1))/ (W(J)-W(J-1)) 
        A = 3.D0* (W(K)-W(K-2)) - (W(J)-W(J-1))* (2.*W(K-1)+W(K+1)) 
        B = -2.0D0* (W(K)-W(K-2)) + (W(J)-W(J-1))* (W(K-1)+W(K+1)) 
        FN(I) = ((B*D+A)*D+W(K-1)* (W(J)-W(J-1)))*D + W(K-2) 
        GN(I) = (3.D0*B*D+2.*A)*D/ (W(J)-W(J-1)) + W(K-1) 
   50 END DO 
   13 J = 2 
      DO 59 I = 1,M 
   60   IF (XD(I)-XN(J)) 61,61,62 
   62   IF (J-N) 76,61,61 
   76   J = J + 1 
        GO TO 60 
   61   D = (XD(I)-XN(J-1))/ (XN(J)-XN(J-1)) 
        A = FN(J-1) + D* ((FN(J-1)-FN(J))*D* (2.0D0*D-3.0D0)+           &
     &      (XN(J)-XN(J-1))* (GN(J-1)* (1.0D0+D* (D-2.0D0))+GN(J)* (D*D-&
     &      D)))                                                        
        RD(I) = YD(I) - A 
   59 END DO 
      IF (IPRINT) 53,53,54 
   54 WRITE (6,FMT=55) 
   55 FORMAT (1H1,24X,39HSPLINE APPROXIMATION OBTAINED BY VB05BD,/,/,4X,&
     &       1HI,10X,5HXN(I),20X,5HFN(I),20X,5HGN(I),/,/)               
      DO 56 I = 1,N 
        WRITE (6,FMT=57) I,XN(I),FN(I),GN(I) 
   57   FORMAT (I5,3E25.14) 
   56 END DO 
      WRITE (6,FMT=58) 
   58 FORMAT (/,/,/,4X,1HI,9X,5HXD(I),18X,5HYD(I),18X,5HWD(I),19X,3HFIT,&
     &       18X,8HRESIDUAL,/,/)                                        
      J = 2 
      DO 159 I = 1,M 
  160   IF (XD(I)-XN(J)) 161,161,162 
  162   IF (J-N) 176,161,161 
  176   J = J + 1 
        WRITE (6,FMT=63) 
   63   FORMAT (5X) 
        GO TO 160 
  161   A = YD(I) - RD(I) 
        WRITE (6,FMT=64) I,XD(I),YD(I),WD(I),A,RD(I) 
   64   FORMAT (I5,5E23.14) 
  159 END DO 
   53 RETURN 
   28 A = 1.0D0/ ((W(I-2)-W(I-1))* (W(I-2)-W(I))* (W(I-2)-W(I+1))*      &
     &    (W(I-2)-W(I+2)))                                              
      B = 1.0D0/ ((W(I-1)-W(I-2))* (W(I-1)-W(I))* (W(I-1)-W(I+1))*      &
     &    (W(I-1)-W(I+2)))                                              
      C = 1.0D0/ ((W(I+1)-W(I-2))* (W(I+1)-W(I-1))* (W(I+1)-W(I))*      &
     &    (W(I+1)-W(I+2)))                                              
      D = 1.0D0/ ((W(I+2)-W(I-2))* (W(I+2)-W(I-1))* (W(I+2)-W(I))*      &
     &    (W(I+2)-W(I+1)))                                              
      GO TO (29,47) IS 
      RETURN 
END                                           
                                                                        
                                                                        
DOUBLE PRECISION FUNCTION TG01BD(IX,N,U,S,D,X) 
 DOUBLE PRECISION X 
 INTEGER IX,N 
 DOUBLE PRECISION D(*),S(*),U(*) 
 DOUBLE PRECISION A,B,EPS,H,Q1,Q2,SS,Z 
 INTEGER IFLG,J 
 DOUBLE PRECISION FD05AD 
 EXTERNAL FD05AD 
 INTRINSIC DABS,DMAX1,MIN0 
 SAVE IFLG,J,H,Q1,Q2,SS,B,A 
 DATA IFLG/0/
      EPS = FD05AD(1)*4.0D0 
      IF (X.LT.U(1)) GO TO 990 
      IF (X.GT.U(N)) GO TO 991 
      IF (IX.LT.0 .OR. IFLG.EQ.0) GO TO 12 
      IF (X.LE.U(J+1)) GO TO 8 
    1 J = J + 1 
   11 IF (X.GT.U(J+1)) GO TO 1 
      GO TO 7 
   12 J = DABS(X-U(1))/ (U(N)-U(1))* (N-1) + 1 
      J = MIN0(J,N-1) 
      IFLG = 1 
      IF (X.GE.U(J)) GO TO 11 
    2 J = J - 1 
      IF (X.LT.U(J)) GO TO 2 
    7 H = U(J+1) - U(J) 
      Q1 = H*D(J) 
      Q2 = H*D(J+1) 
      SS = S(J+1) - S(J) 
      B = 3D0*SS - 2D0*Q1 - Q2 
      A = Q1 + Q2 - 2D0*SS 
    8 Z = (X-U(J))/H 
      TG01BD = ((A*Z+B)*Z+Q1)*Z + S(J) 
      RETURN 
  990 IF (X.LE.U(1)-EPS*DMAX1(DABS(U(1)),DABS(U(N)))) GO TO 99 
      J = 1 
      GO TO 7 
  991 IF (X.GE.U(N)+EPS*DMAX1(DABS(U(1)),DABS(U(N)))) GO TO 99 
      J = N - 1 
      GO TO 7 
   99 IFLG = 0 
      TG01BD = 0D0 
      RETURN 
END                                           
                                                                        
                                                                        
DOUBLE PRECISION FUNCTION FD05AD(INUM) 
 INTEGER INUM 
 DOUBLE PRECISION DC(5) 
 SAVE DC 
 DATA DC(1)/2.2204460492504D-16/ 
 DATA DC(2)/1.1102230246253D-16/ 
 DATA DC(4)/2.2250738585073D-308/ 
 DATA DC(5)/1.7976931348622D+308/ 
                                                                        
      IF ( INUM .LE. 0 ) THEN 
         FD05AD = DC( 1 ) 
      ELSE IF ( INUM .GE. 6 ) THEN 
         FD05AD = DC( 5 ) 
      ELSE IF ( INUM .EQ. 3 ) THEN 
         FD05AD = DC(4)/2.0D0**52 
      ELSE 
         FD05AD = DC( INUM ) 
      ENDIF 
      RETURN 
END                                
