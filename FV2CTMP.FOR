// JOB 'NZZ123,123456789'
//  EXEC CRAY
//SYSUT1  DD  *
# USER=osu310 PW=hccsc3
#@$-r ts5964c
# @$-eo
# @$-lT 300
# @$-lM 1.5mw
# @$
set -x
cat << 'EOF' > for.f
      PROGRAM FLOW
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                              C
C     *****FLOWV2.FOR*****      15/04/88   VERSION 2.0                         C
C                                                                              C
C     TWO DIMEMSIONAL ICE FLOW MODEL                                           C
C     LAMELLAR FLOW                                                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION H(200,200),U(200,200),V(200,200),
     1UNEW(200,200),VNEW(200,200),QNEWX(200,200),QNEWY(200,200),
     2DELTHK(200,200),OLDH(200,200),D(200,200),
     3BOTTOM(200)
      PARAMETER (BDOT=9.6D-09,A=5.2D-25,TB=0.1D06,XN=3.0D0,
     1T=3.1536D07,RHO=1.D03,G=9.8D0,DXL=1.0D06,SECYR=3.1536D07,
     2DYL=1.0D06,NX=10,NY=10)
CCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC
C ZERO ALL VARIABLES C
CCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC
      ITEST=0
      X=0.0D0
      LIM=0
      SAVE=0.0D0
      CHGMAX=0.0D0
      DO 8 I=1,101
      DO 9 J=1,101
      DELTHK(I,J)=0.0D0
      U(I,J)=0.0D0
      V(I,J)=0.0D0
      H(I,J)=0.0D0
      UNEW(I,J)=0.0D0
      VNEW(I,J)=0.0D0
      QNEWX(I,J)=0.0D0
      OLDH(I,J)=0.0D0
      D(I,J)=0.0D0
9     CONTINUE
8     CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCC
C SET INITIAL CONDITIONS C
CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCC
      DT=T/4.0D0
      TOL=2.0D0
      DX=DXL/FLOAT(NX)
      DY=DYL/FLOAT(NY)
      DXNEW=1.0D0/(2.0D0*DX)
      DYNEW=1.0D0/(2.0D0*DY)
      XN2=1.0D0/(XN+2.0D0)
      DY2=1.0D0/(2.0D0*DY**2)
      DX2=1.0D0/(2.0D0*DX**2)
      DX2A=1.0D0/(DX**2)
      RG=1.0D0/(RHO*G)
C-------SET UP INITIAL PROFILE---------------------------------------C
      CALL INIT(H,OLDH,RG,TB,X,DX,DXL,NX,NY,BOTTOM)
CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCC
C PERFORM ITERATION      C
CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCC
      CALL VEL(U,V,NX,NY,H,A,RHO,G,DXNEW,DYNEW,XN,XN2)
C---------ITERATION STATEMENT----------------------------------------C
50    IF(TOL.LE.1.0D-05.OR.LIM.GE.40000.OR.ITEST.NE.0) GOTO 100
      TOL=0.0D0
      LIM=LIM+1
      DO 1 I=1,NX+1
      DO 2 J=1,NY+1
      OLDH(I,J)=H(I,J)
2     CONTINUE
1     CONTINUE
      CALL DIFFCOE(NX,NY,H,A,RHO,G,XN,XN2,DXNEW,DYNEW,D)
      CALL HIEGHT(H,OLDH,D,DT,DX2,DX2A,DY2,BDOT,BOTTOM,NX,NY)
      CALL VEL(U,V,NX,NY,H,A,RHO,G,DXNEW,DYNEW,XN,XN2)
      DO 17 I=1,NX+1
      DO 13 J=1,NY+1
      DELTHK(I,J)=H(I,J)-OLDH(I,J)
      CHGMAX=DABS(DELTHK(I,J))
      IF(TOL.LT.CHGMAX) TOL=CHGMAX
13    CONTINUE
17    CONTINUE
C---------TEST SECTION-----------------------------------------------C
      DO 21 I=1,NX+1
      DO 25 J=1,NY+1
      UNEW(I,J)=U(I,J)*SECYR
      VNEW(I,J)=V(I,J)*SECYR
      QNEWX(I,J)=U(I,J)*H(I,J)*(1.0D0)*SECYR
      QNEWY(I,J)=V(I,J)*H(I,J)*(1.0D0)*SECYR
C      IF(UNEW(I,J).LT.0.0D0.OR.UNEW(I,J).GT.2.0D03) ITEST=1
C      IF(I.NE.NX+1.AND.(H(I,J)-H(I+1,J)).LT.0.0D0) ITEST=2
25    CONTINUE
21    CONTINUE
      GOTO 50
100   CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCC
C WRITE THE OUTPUT       C
CCCCCCCCCCCCCCCCCCCCCCCCCC
      NNY=NY/10
      WRITE(6,903)
903   FORMAT('\n','THE NEW ICE THICKNESS')
      DO 14 I=1,NX+1
      WRITE(6,901) (H(I,J),J=1,NY+1,NNY)
901   FORMAT('\n',11(1X,F9.4))
14    CONTINUE
      WRITE(6,913)
913   FORMAT('\n','THE OLDH ICE THICKNESS')
      DO 84 I=1,NX+1
      WRITE(6,901) (OLDH(I,J),J=1,NY+1,NNY)
84    CONTINUE
      WRITE(6,906)
906   FORMAT('\n','THE CHANGE IN THICKNESS')
      DO 92 I=1,NX+1
      WRITE(6,901) (DELTHK(I,J),J=1,NY+1,NNY)
92    CONTINUE
      WRITE(6,905)
905   FORMAT('\n\n\n','THE  U  VELOCITIES')
      DO 27 I=1,NX+1
      WRITE(6,901) (UNEW(I,J),J=1,NY+1,NNY)
27    CONTINUE
C      WRITE(6,907)
C907   FORMAT('\n\n\n','THE  V  VELOCITIES')
C      DO 19 I=1,NX+1
C      WRITE(6,901) (VNEW(I,J),J=1,NY+1,NNY)
C19    CONTINUE
      WRITE(6,912)
912   FORMAT('\n','THE DIFFUSIVITY COEFFICIENT')
      DO 81 I=1,NX+1
      WRITE(6,910) (D(I,J),J=1,NY+1,NNY)
81    CONTINUE
      WRITE(6,908)
908   FORMAT('\n\n\n','THE DISCHARGE IN THE X-DIRECTION')
      DO 79 I=1,NX+1
      WRITE(6,910) (QNEWX(I,J),J=1,NY+1,NNY)
910   FORMAT('\n',11(1X,D11.6))
79    CONTINUE
C      WRITE(6,909)
C909   FORMAT('\n\n\n','THE DISCHARGE IN THE Y-DIRECTION')
C      DO 80 I=1,NX+1
C      WRITE(6,910) (QNEWY(I,J),J=1,NY+1,NNY)
C80    CONTINUE
      WRITE(6,902) LIM
902   FORMAT('\n','THE NUMBER OF ITERATIONS =',I8)
      WRITE(6,904) ITEST
904   FORMAT('\n','ITEST=',I3)
      STOP
      END
      SUBROUTINE VEL(U,V,NX,NY,H,A,RHO,G,DXNEW,DYNEW,XN,XN2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DIMENSION U(200,200),V(200,200),H(200,200)
      DO 200 I=2,NX
      DO 201 J=2,NY
      U(I,J)=2.0D0*A*(RHO*G)**XN*(((H(I+1,J)-H(I-1,J))*DXNEW)**2+((
     1H(I,J+1)-H(I,J-1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,J)**(
     2XN+1.0D0))*XN2*DABS(H(I+1,J)-H(I-1,J))*DXNEW
      V(I,J)=2.0D0*A*(RHO*G)**XN*(((H(I+1,J)-H(I-1,J))*DXNEW)**2+((
     1H(I,J+1)-H(I,J-1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,J)**(
     2XN+1.0D0))*XN2*(H(I,J+1)-H(I,J-1))*DYNEW
201   CONTINUE
200   CONTINUE
C*********ICE DIVIDE, END OF GLACIER LOOPS***************************C
      DO 202 J=2,NY
      U(1,J)=0.0D0
      V(1,J)=0.0D0
      U(NX+1,J)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,J)-4.0D0*H(NX,J)+3.0D0*H(N
     1X+1,J))*DXNEW)**2+((H(NX+1,J+1)-H(NX+1,NY+1))*DYNEW)**2)**((XN-1.0
     2D0)*0.5D0)*(H(NX+1,J)**(XN+1.0D0))*XN2*DABS(H(NX-1,J)-4.0D0*H(NX
     3,J)+3.0D0*H(NX+1,J))*DXNEW
      V(NX+1,J)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,J)-4.0D0*H(NX,J)+3.0D0*H(N
     1X+1,J))*DXNEW)**2+((H(NX+1,J+1)-H(NX+1,NY+1))*DYNEW)**2)**((XN-1.0
     2D0)*0.5D0)*H(NX+1,J)**(XN+1.0D0)*XN2*(H(NX+1,J+1)-H(NX+1,J-1))*DYN
     3EW
202   CONTINUE
C**********TOP,BOTTOM LOOPS******************************************C
      DO 203 I=2,NX
      U(I,1)=2.0D0*A*(RHO*G)**XN*(((H(I+1,1)-H(I-1,1))*DXNEW)**2+((
     1H(I,2)-H(I,NY+1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,1)**(
     2XN+1.0D0))*XN2*DABS(H(I+1,1)-H(I-1,1))*DXNEW
      V(I,1)=2.0D0*A*(RHO*G)**XN*(((H(I+1,1)-H(I-1,1))*DXNEW)**2+((
     1H(I,2)-H(I,NY+1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,1)**(XN+1
     2.0D0))*XN2*(H(I,2)-H(I,NY+1))*DYNEW
      U(I,NY+1)=2.0D0*A*(RHO*G)**XN*(((H(I+1,NY+1)-H(I-1,NY+1))*DXNEW)**
     12+((H(I,1)-H(I,NY))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,NY+1)**(XN
     2+1.0D0))*XN2*DABS(H(I+1,NY+1)-H(I-1,NY+1))*DXNEW
      V(I,NY+1)=2.0D0*A*(RHO*G)**XN*(((H(I+1,NY+1)-H(I-1,NY+1))*DXNEW)**
     12+((H(I,1)-H(I,NY))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,NY+1)**(XN
     2+1.0D0))*XN2*(H(I,1)-H(I,NY))*DYNEW
203   CONTINUE
C**********THE 4 CORNERS*********************************************C
      U(1,1)=0.0D0
      V(1,1)=0.0D0
      U(NX+1,1)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,1)-4.0D0*H(NX,1)+3.0D0*H(N
     1X+1,1))*DXNEW)**2+((H(NX+1,2)-H(NX+1,NY+1))*DYNEW)**2)**((XN-1.0D
     20)*0.5D0)*(H(NX+1,1)**(XN+1.0D0))*XN2*DABS(H(NX-1,1)-4.0D0*H(NX,
     31)+3.0D0*H(NX+1,1))*DXNEW
      V(NX+1,1)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,1)-4.0D0*H(NX,1)+3.0D0*H(N
     1X+1,1))*DXNEW)**2+((H(NX+1,2)-H(NX+1,NY+1))*DYNEW)**2)**((XN-1.0D
     20)*0.5D0)*H(NX+1,1)**(XN+1.0D0)*XN2*(H(NX+1,2)-H(NX+1,NY+1))*DYNEW
      U(1,NY+1)=0.0D0
      V(1,NY+1)=0.0D0
      U(NX+1,NY+1)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,NY+1)-4.0D0*H(NX,NY+1)+
     13.0D0*H(NX+1,NY+1))*DXNEW)**2+((H(NX+1,1)-H(NX+1,NY))*DYNEW)**2)**
     2((XN-1.0D0)*0.5D0)*(H(NX+1,NY+1)**(XN+1.0D0))*XN2*DABS(H(NX-1,NY+1
     3)-4.0D0*H(NX,NY+1)+3.0D0*H(NX+1,NY+1))*DXNEW
      V(NX+1,NY+1)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,NY+1)-4.0D0*H(NX,NY+1)+
     13.0D0*H(NX+1,NY+1))*DXNEW)**2+((H(NX+1,1)-H(NX+1,NY))*DYNEW)**2)**
     2((XN-1.0D0)*0.5D0)*H(NX+1,NY+1)**(XN+1.0D0)*XN2*(H(NX+1,1)-H(NX+1,
     3NY))*DYNEW
      RETURN
      END
      SUBROUTINE DIFFCOE(NX,NY,H,A,RHO,G,XN,XN2,DXNEW,DYNEW,D)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DIMENSION H(200,200),D(200,200)
      DO 200 I=2,NX
      DO 201 J=2,NY
      D(I,J)=(2.0D0*A*(RHO*G)**(XN)*(H(I,J)**(XN+2.0D0))*
     1(((H(I+1,J)-H(I-1,J))*DXNEW)**2+((H(I,J+1)-H(I,J-1))*DYNEW
     2)**2)**((XN-1.0D0)*0.5D0))*XN2
201   CONTINUE
200   CONTINUE
C*********ICE DIVIDE, END OF GLACIER LOOPS***************************C
      DO 202 J=2,NY
      D(1,J)=2.0D0*A*(RHO*G)**XN*H(1,J)**(XN+2.0D0)*(((H(1,J+1)-H(1,J-1
     1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*XN2
      D(NX+1,J)=2.0D0*A*(RHO*G)**XN*H(NX+1,J)**(XN+2.0D0)*(((H(NX-1,J)-4
     1.0D0*H(NX,J)+3.0D0*H(NX+1,J))*DXNEW)**2+((H(NX+1,J+1)-H(NX+1,J-1))
     2*DYNEW)**2)**((XN-1.0D0)*0.5D0)*XN2
202   CONTINUE
C**********TOP,BOTTOM LOOPS******************************************C
      DO 203 I=2,NX
      D(I,1)=(2.0D0*A*(RHO*G)**(XN)*(H(I,1)**(XN+2.0D0))*
     1(((H(I+1,1)-H(I-1,1))*DXNEW)**2+((H(I,2)-H(I,NY+1))*DYNEW
     2)**2)**((XN-1.0D0)*0.5D0))*XN2
      D(I,NY+1)=(2.0D0*A*(RHO*G)**(XN)*(H(I,NY+1)**(XN+2.0D0))*
     1(((H(I+1,NY+1)-H(I-1,NY+1))*DXNEW)**2+((H(I,1)-H(I,NY))*DYNEW)
     2**2)**((XN-1.0D0)*0.5D0))*XN2
203   CONTINUE
C**********THE 4 CORNERS*********************************************C
      D(1,1)=2.0D0*A*(RHO*G)**XN*H(1,1)**(XN+2.0D0)*(((H(1,2)-H(1,NY+1))
     1*DYNEW)**2)**((XN-1.0D0)*0.5D0)*XN2
      D(NX+1,1)=2.0D0*A*(RHO*G)**XN*H(NX+1,1)**(XN+2.0D0)*(((H(NX-1,1)-4
     1.0D0*H(NX,1)+3.0D0*H(NX+1,1))*DXNEW)**2+((H(NX+1,2)-H(NX+1,NY+1))*
     2DYNEW)**2)**((XN-1.0D0)*0.5D0)*XN2
      D(1,NY+1)=2.0D0*A*(RHO*G)**XN*H(1,NY+1)**(XN+2.0D0)*(((H(1,1)-H(1,
     1NY))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*XN2
      D(NX+1,NY+1)=2.0D0*A*(RHO*G)**XN*H(NX+1,NY+1)**(XN+2.0D0)*(((H(NX-
     11,NY+1)-4.0D0*H(NX,NY+1)+3.0D0*H(NX+1,NY+1))*DXNEW)**2+((H(NX+1,1)
     2-H(NX+1,NY))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*XN2
      RETURN
      END
      SUBROUTINE INIT(H,OLDH,RG,TB,X,DX,DXL,NX,NY,BOTTOM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DIMENSION H(200,200),OLDH(200,200),BOTTOM(200)
      XXX=2.0D0*DX/10.0D0
      DO 101 J=1,NY+1
      X=0.0D0
      DO 10 I=1,NX+1
      H(I,J)=DSQRT(2.0D0*TB*RG*(DXL+XXX-X))
      OLDH(I,J)=H(I,J)
      X=X+DX
10    CONTINUE
101   CONTINUE
      DO 91 J=1,NY+1
      BOTTOM(J)=H(NX+1,J)
91    CONTINUE
      RETURN
      END
      SUBROUTINE HIEGHT(H,OLDH,D,DT,DX2,DX2A,DY2,BDOT,BOTTOM,NX,NY)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DIMENSION OLDH(200,200),H(200,200),D(200,200),BOTTOM(200)
C*********MAIN BODY LOOPS********************************************C
      DO 300 I=2,NX
      DO 301 J=2,NY
      H(I,J)=(((D(I+1,J)+D(I,J))*(OLDH(I+1,J)-OLDH(I,J))
     1-(D(I,J)+D(I-1,J))*(OLDH(I,J)-OLDH(I-1,J)))*DX2
     2+((D(I,J+1)+D(I,J))*(OLDH(I,J+1)-OLDH(I,J))-(D(I,J)+D(I,J-
     31))*(OLDH(I,J)-OLDH(I,J-1)))*DY2+BDOT)*DT+OLDH(I,J)
301   CONTINUE
300   CONTINUE
C**********ICE DIVIDE, END OF GLACIER LOOP***************************C
      DO 302 J=2,NY
      H(1,J)=((D(2,J)+D(1,J))*(OLDH(2,J)-OLDH(1,J))*DX2A+
     1(D(1,J+1)+D(1,J))*(OLDH(1,J+1)-OLDH(1,J))-(D(1,J)+D(1,J-1)
     2)*(OLDH(1,J)-OLDH(1,J-1))*DY2+BDOT)*DT+OLDH(1,J)
      H(NX+1,J)=BOTTOM(J)
302   CONTINUE
C**********TOP, BOTTOM LOOPS*****************************************C
      DO 303 I=2,NX
      H(I,1)=(((D(I+1,1)+D(I,1))*(OLDH(I+1,1)-OLDH(I,1))
     1-(D(I,1)+D(I-1,1))*(OLDH(I,1)-OLDH(I-1,1)))*DX2
     2+((D(I,2)+D(I,1))*(OLDH(I,2)-OLDH(I,1))-(D(I,1)+D(I,NY+1))
     3*(OLDH(I,1)-OLDH(I,NY+1)))*DY2+BDOT)*DT+OLDH(I,1)
      H(I,NY+1)=(((D(I+1,NY+1)+D(I,NY+1))*(OLDH(I+1,NY+1)-OLDH(I,NY+
     11))-(D(I,NY+1)+D(I-1,NY+1))*(OLDH(I,NY+1)-OLDH(I-1,NY+1)))
     2*DX2
     3+((D(I,1)+D(I,NY+1))*(OLDH(I,1)-OLDH(I,NY+1))-(D(I,NY+1)
     4+D(I,NY))
     5*(OLDH(I,NY+1)-OLDH(I,NY)))*DY2+BDOT)*DT+OLDH(I,NY+1)
303   CONTINUE
C**********THE 4 CORNERS*********************************************C
      H(1,1)=((D(2,1)+D(1,1))*(OLDH(2,1)-OLDH(1,1))*DX2A+
     1((D(I,2)+D(1,1))*(OLDH(1,2)-OLDH(1,1))-(D(1,1)+D(1,NY+1))*
     2(OLDH(1,1)-OLDH(1,NY+1)))*DY2+BDOT)*DT+OLDH(1,1)
      H(1,NY+1)=((D(2,NY+1)+D(1,NY+1))*(OLDH(2,NY+1)-OLDH(1,NY+1))
     1*DX2A+
     2((D(1,1)+D(1,NY+1))*(OLDH(1,1)-OLDH(1,NY+1))-(D(1,NY+1)
     3+D(1,NY))
     4*(OLDH(1,NY+1)-OLDH(1,NY)))*DY2+BDOT)*DT+OLDH(1,NY+1)
      H(NX+1,1)=BOTTOM(J)
      H(NX+1,NY+1)=BOTTOM(J)
      RETURN
      END
      SUBROUTINE CHK(NX,NY)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(6,200) NX,NY
200   FORMAT('\n','NX= ',I4,5X,'NY= ',I4)
      RETURN
      END
EOF
time cft -dp -eL for.f
cat for.l
segldr for.o
time a.out
balance
rm for.f for.l for.o a.out
//SYSIN  DD  *
CRSUBMIT FILE(SYSUT1) COLS(1 72) DEBLANK HOLD
