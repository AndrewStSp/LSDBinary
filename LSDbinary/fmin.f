! Ax,Bx - min&max interval for x
! F -external function f(x)
! Tol - tolerantnost (1.e-5)
      real(8) FUNCTION FMIN(AX,BX,F,TOL,Fout,number)
	implicit real(8) (a-h,o-z)
	real Fout
      C=0.5d0*(3.d0-SQRT(5.d0))
      EPS=1.d0
10    EPS=EPS/2.d0
      TOL1=1.d0+EPS
      IF(TOL1.GT.1.d0)GO TO 10
      EPS=SQRT(EPS)
      A=AX
      B=BX
      V=A+C*(B-A)
      W=V
      X=V
      E=0.
     	FX=F(X,number)
     
	FV=FX
      FW=FX
20     XM=.5d0*(A+B)
      TOL1=EPS*ABS(X)+TOL/3.d0
      TOL2=2.*TOL1
       IF(ABS(X-XM).LE.(TOL2-.5d0*(B-A)))GO TO 90
      IF(ABS(E).LE.TOL1)GO TO 40
      R=(X-W)*(FX-FV)
      Q=(X-V)*(FX-FW)
      P=(X-V)*Q-(X-W)*R
      Q=2.*(Q-R)
      IF(Q.GT.0.)P=-P
      Q=ABS(Q)
      R=E
      E=D
30    IF(ABS(P).GE.ABS(.5d0*Q*R))GO TO 40
      IF(P.LE.Q*(A-X))GO TO 40
      IF(P.GE.Q*(B-X))GO TO 40
      D=P/Q
      U=X+D
       IF((U-A).LT.TOL2)D=SIGN(TOL1,XM-X)
       IF((B-U).LT.TOL2)D=SIGN(TOL1,XM-X)
      GO TO 50
40    IF(X.GE.XM)E=A-X
      IF(X.LT.XM)E=B-X
      D=C*E
50    IF(ABS(D).GE.TOL1)U=X+D
      IF(ABS(D).LT.TOL1)U=X+SIGN(TOL1,D)
      FU=F(U,number)
      IF(FU.GT.FX)GO TO 60
      IF(U.GE.X)A=X
      IF(U.LT.X)B=X
       V=W
      FV=FW
      W=X
      FW=FX
      X=U
      FX=FU
      GO TO 20
60    IF(U.LT.X)A=U
       IF(U.GE.X)B=U
       IF(FU.LE.FW)GO TO 70
      IF(W.EQ.X)GO TO 70
       IF(FU.LE.FV)GO TO 80
      IF(V.EQ.X)GO TO 80
        IF(V.EQ.W)GO TO 80
      GO TO 20
70    V=W
      FV=FW
      W=U
      FW=FU
      GO TO 20
80    V=U
      FV=FU
      GO TO 20
90    FMIN=X
	Fout=FU
      RETURN
      END
