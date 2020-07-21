// JOB ,
// TIME=(5,00),REGION=1024K
/*JOBPARM LINES=3000
//  EXEC VSF2CLG,TIME=(5,00)
//FORT.SYSIN  DD  *
      PROGRAM FLOW
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                              C
C     *****FLOWV2N1.FOR*****      16/10/88   VERSION 1.0                       C
C                                                                              C
C     TWO DIMEMSIONAL ICE FLOW MODEL                                           C
C     USING OERLEMANS DIFFUSION SCHEME                                         C
C     MODIFIED FOR THE CRAY AND IN SUBROUTINE FORM                             C
C     LAMELLAR AND NORMAL FLOW                                                 C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
C     DIMENSIONED COMMON BLOCKS
      COMMON /A/ H(200,200),U(200,200),V(200,200),D(200,200)
      COMMON /B/ OLDH(200,200),UNEW(200,200),VNEW(200,200)
      COMMON /C/ QNEWX(200,200),QNEWY(200,200),DELTHK(200,200),
     1           BOTTOM(200)
      COMMON /C3/ DXX(200,200),DXY(200,200)
      COMMON /C4/ OLDU(200,200),OLDV(200,200),OLDDXX(200,200),OLDDXY(200
     1,200),UB(200,200),VB(200,200)
C     PARAMETER COMMON BLOCKS
      COMMON /D/ BDOT,A,TB,XN,DT,RHO,G,DXL,DYL,SECYR
C     CALCULATED COMMON BLOCKS
      COMMON /E/ NX,NY,DXNEW,DYNEW,XN2,RG,X,DX,DY
C     PARAMETER (BDOT=9.6D-09,A=5.2D-25,TB=0.1D06,XN=3.0D0,
C    1T=3.1536D07,RHO=1.D03,G=9.8D0,DXL=1.0D06,SECYR=3.1536D07,
C    2DYL=1.0D06,NX=10,NY=10)
CCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC
C ZERO ALL VARIABLES C
CCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCC
      BDOT=9.6D-09
      A=5.2D-25
      TB=0.1D06
      XN=3.0D0
      T=3.1536D07
      RHO=1.0D03
      G=9.8D0
      DXL=1.0D06
      SECYR=3.1536D07
      DYL=1.0D06
      NX=10
      NY=10
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
      DXX(I,J)=0.0D0
      DXY(I,J)=0.0D0
      D(I,J)=0.0D0
      UB(I,J)=0.0D0
      VB(I,J)=0.0D0
9     CONTINUE
8     CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCC
C SET INITIAL CONDITIONS C
CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCC
C      DT=T/4.0D0
      DT=T
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
C-------SET UP INITIAL PROFILE AND VELOCITIES------------------------C
      CALL INIT
      CALL VELINT
      CALL SIGINT
CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCC
C PERFORM ITERATION      C
CCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCC
C---------ITERATION STATEMENT----------------------------------------C
50    IF(TOL.LE.1.0D-05.OR.LIM.GE.0.OR.ITEST.NE.0) GOTO 100
      TOL=0.0D0
      LIM=LIM+1
      DO 1 I=1,NX+1
      DO 2 J=1,NY+1
      OLDH(I,J)=H(I,J)
2     CONTINUE
1     CONTINUE
      CALL VEL
      CALL STRESS
      CALL VEL
      CALL STRESS
C      COMPUTE DIFFUSION COEFFICIENT                                   C
      CALL DIFFCOE
C      COMPUTE ICE SHEET ELEVATIONS                                    C
      CALL HIEGHT
      CALL VEL
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
      IF(UNEW(I,J).LT.0.0D0.OR.UNEW(I,J).GT.2.0D03) ITEST=1
      IF(I.NE.NX+1.AND.(H(I,J)-H(I+1,J)).LT.0.0D0) ITEST=2
25    CONTINUE
21    CONTINUE
      GOTO 50
100   CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCC
C WRITE THE OUTPUT       C
CCCCCCCCCCCCCCCCCCCCCCCCCC
      NNY=NY/10
      WRITE(6,903)
903   FORMAT(' ','THE NEW ICE THICKNESS')
      DO 14 I=1,NX+1
      WRITE(6,901) (H(I,J),J=1,NY+1,NNY)
901   FORMAT(' ',11(1X,F9.4))
14    CONTINUE
      WRITE(6,913)
913   FORMAT(' ','THE OLDH ICE THICKNESS')
      DO 84 I=1,NX+1
      WRITE(6,901) (OLDH(I,J),J=1,NY+1,NNY)
84    CONTINUE
      WRITE(6,906)
906   FORMAT(' ','THE CHANGE IN THICKNESS')
      DO 92 I=1,NX+1
      WRITE(6,901) (DELTHK(I,J),J=1,NY+1,NNY)
92    CONTINUE
      WRITE(6,905)
905   FORMAT(' ','THE  U  VELOCITIES')
      DO 27 I=1,NX+1
      WRITE(6,901) (UNEW(I,J),J=1,NY+1,NNY)
27    CONTINUE
C      WRITE(6,907)
C907   FORMAT(' ','THE  V  VELOCITIES')
C      DO 19 I=1,NX+1
C      WRITE(6,901) (VNEW(I,J),J=1,NY+1,NNY)
C19    CONTINUE
      WRITE(6,912)
912   FORMAT(' ','THE DIFFUSIVITY COEFFICIENT')
      DO 81 I=1,NX+1
      WRITE(6,910) (D(I,J),J=1,NY+1,NNY)
81    CONTINUE
      WRITE(6,908)
908   FORMAT(' ','THE DISCHARGE IN THE X-DIRECTION')
      DO 79 I=1,NX+1
      WRITE(6,910) (QNEWX(I,J),J=1,NY+1,NNY)
910   FORMAT(' ',11(1X,D11.6))
79    CONTINUE
C      WRITE(6,909)
C909   FORMAT(' ','THE DISCHARGE IN THE Y-DIRECTION')
C      DO 80 I=1,NX+1
C      WRITE(6,910) (QNEWY(I,J),J=1,NY+1,NNY)
C80    CONTINUE
      WRITE(6,902) LIM
902   FORMAT(' ','THE NUMBER OF ITERATIONS =',I8)
      WRITE(6,904) ITEST
904   FORMAT(' ','ITEST=',I3)
      STOP
      END
      SUBROUTINE VEL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CHGU(200,200),CHGV(200,200)
C     DIMENSIONED COMMON BLOCKS
      COMMON /A/ H(200,200),U(200,200),V(200,200),D(200,200)
      COMMON /B/ OLDH(200,200),UNEW(200,200),VNEW(200,200)
      COMMON /C/ QNEWX(200,200),QNEWY(200,200),DELTHK(200,200),
     1           BOTTOM(200)
      COMMON /C3/ DXX(200,200),DXY(200,200)
      COMMON /C4/ OLDU(200,200),OLDV(200,200),OLDDXX(200,200),OLDDXY(200
     1,200),UB(200,200),VB(200,200)
C     PARAMETER COMMON BLOCKS
      COMMON /D/ BDOT,A,TB,XN,DT,RHO,G,DXL,DYL,SECYR
C     CALCULATED COMMON BLOCKS
      COMMON /E/ NX,NY,DXNEW,DYNEW,XN2,RG,X,DX,DY
      DO 300 I=1,NX+1
      DO 301 J=1,NY+1
      OLDU(I,J)=U(I,J)
      OLDV(I,J)=V(I,J)
301   CONTINUE
300   CONTINUE
C******MAIN BODY LOOPS*************************************************C
      DO 200 I=2,NX
      DO 201 J=2,NY
      U(I,J)=2.0D0*A*(RHO*G)**XN*(((H(I+1,J)-H(I-1,J))*DXNEW)**2+((
     1H(I,J+1)-H(I,J-1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,J)**(
     2XN+1.0D0))*XN2*DABS(H(I+1,J)-H(I-1,J))*DXNEW
     3+2.0D0*A*(DXX(I,J)**2+DXY(I,J)**2)*(-RHO*G)*H(I,J)**2*(H(I+1,J)-H(
     4I-1,J))*DXNEW/XN+UB(I,J)
      V(I,J)=2.0D0*A*(RHO*G)**XN*(((H(I+1,J)-H(I-1,J))*DXNEW)**2+((
     1H(I,J+1)-H(I,J-1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,J)**(
     2XN+1.0D0))*XN2*(H(I,J+1)-H(I,J-1))*DYNEW
     3+2.0D0*A*(DXX(I,J)**2+DXY(I,J)**2)*(-RHO*G)*H(I,J)**2*(H(I,J+1)-H(
     4I,J-1))*DYNEW/XN+VB(I,J)
201   CONTINUE
200   CONTINUE
C*********ICE DIVIDE, END OF GLACIER LOOPS***************************C
      DO 202 J=2,NY
      U(1,J)=0.0D0+UB(1,J)
      V(1,J)=0.0D0+VB(1,J)
      U(NX+1,J)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,J)-4.0D0*H(NX,J)+3.0D0*H(N
     1X+1,J))*DXNEW)**2+((H(NX+1,J+1)-H(NX+1,NY+1))*DYNEW)**2)**((XN-1.0
     2D0)*0.5D0)*(H(NX+1,J)**(XN+1.0D0))*XN2*DABS(H(NX-1,J)-4.0D0*H(NX
     3,J)+3.0D0*H(NX+1,J))*DXNEW
     3+2.0D0*A*(DXX(NX+1,J)**2+DXY(NX+1,J)**2)*(-RHO*G)*H(NX+1,J)**2*(H(
     4NX-1,J)-4.0D0*H(NX,J)+3.0D0*H(NX+1,J))*DXNEW/XN+UB(NX+1,J)
      V(NX+1,J)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,J)-4.0D0*H(NX,J)+3.0D0*H(N
     1X+1,J))*DXNEW)**2+((H(NX+1,J+1)-H(NX+1,NY+1))*DYNEW)**2)**((XN-1.0
     2D0)*0.5D0)*H(NX+1,J)**(XN+1.0D0)*XN2*(H(NX+1,J+1)-H(NX+1,J-1))*DYN
     3EW
     3+2.0D0*A*(DXX(NX+1,J)**2+DXY(NX+1,J)**2)*(-RHO*G)*H(NX+1,J)**2*(H(
     4NX+1,J+1)-H(NX+1,J-1))*DYNEW/XN+VB(NX+1,J)
202   CONTINUE
C**********TOP,BOTTOM LOOPS******************************************C
      DO 203 I=2,NX
      U(I,1)=2.0D0*A*(RHO*G)**XN*(((H(I+1,1)-H(I-1,1))*DXNEW)**2+((
     1H(I,2)-H(I,NY+1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,1)**(
     2XN+1.0D0))*XN2*DABS(H(I+1,1)-H(I-1,1))*DXNEW
     3+2.0D0*A*(DXX(I,1)**2+DXY(I,1)**2)*(-RHO*G)*H(I,1)**2*(H(I+1,1)-H(
     4I-1,1))*DXNEW/XN+UB(I,1)
      V(I,1)=2.0D0*A*(RHO*G)**XN*(((H(I+1,1)-H(I-1,1))*DXNEW)**2+((
     1H(I,2)-H(I,NY+1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,1)**(XN+1
     2.0D0))*XN2*(H(I,2)-H(I,NY+1))*DYNEW
     3+2.0D0*A*(DXX(I,1)**2+DXY(I,1)**2)*(-RHO*G)*H(I,1)**2*(H(I,2)-H(
     4I,NY+1))*DYNEW/XN+VB(I,1)
      U(I,NY+1)=2.0D0*A*(RHO*G)**XN*(((H(I+1,NY+1)-H(I-1,NY+1))*DXNEW)**
     12+((H(I,1)-H(I,NY))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,NY+1)**(XN
     2+1.0D0))*XN2*DABS(H(I+1,NY+1)-H(I-1,NY+1))*DXNEW
     3+2.0D0*A*(DXX(I,NY+1)**2+DXY(I,NY+1)**2)*(-RHO*G)*H(I,NY+1)**2*(H(
     4I+1,NY+1)-H(I-1,NY+1))*DXNEW/XN+UB(I,NY+1)
      V(I,NY+1)=2.0D0*A*(RHO*G)**XN*(((H(I+1,NY+1)-H(I-1,NY+1))*DXNEW)**
     12+((H(I,1)-H(I,NY))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,NY+1)**(XN
     2+1.0D0))*XN2*(H(I,1)-H(I,NY))*DYNEW
     3+2.0D0*A*(DXX(I,NY+1)**2+DXY(I,NY+1)**2)*(-RHO*G)*H(I,NY+1)**2*(H(
     4I,1)-H(I,NY))*DYNEW/XN+VB(I,NY+1)
203   CONTINUE
C**********THE 4 CORNERS*********************************************C
      U(1,1)=0.0D0+UB(1,1)
      V(1,1)=0.0D0+VB(1,1)
      U(NX+1,1)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,1)-4.0D0*H(NX,1)+3.0D0*H(N
     1X+1,1))*DXNEW)**2+((H(NX+1,2)-H(NX+1,NY+1))*DYNEW)**2)**((XN-1.0D
     20)*0.5D0)*(H(NX+1,1)**(XN+1.0D0))*XN2*DABS(H(NX-1,1)-4.0D0*H(NX,
     31)+3.0D0*H(NX+1,1))*DXNEW
     3+2.0D0*A*(DXX(NX+1,1)**2+DXY(NX+1,1)**2)*(-RHO*G)*H(NX+1,1)**2*(H(
     4NX-1,1)-4.0D0*H(NX,1)+3.0D0*H(NX+1,1))*DXNEW/XN+UB(NX+1,1)
      V(NX+1,1)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,1)-4.0D0*H(NX,1)+3.0D0*H(N
     1X+1,1))*DXNEW)**2+((H(NX+1,2)-H(NX+1,NY+1))*DYNEW)**2)**((XN-1.0D
     20)*0.5D0)*H(NX+1,1)**(XN+1.0D0)*XN2*(H(NX+1,2)-H(NX+1,NY+1))*DYNEW
     3+2.0D0*A*(DXX(NX+1,1)**2+DXY(NX+1,1)**2)*(-RHO*G)*H(NX+1,1)**2*(H(
     4NX+1,2)-H(NX+1,NY+1))*DYNEW/XN+VB(NX+1,1)
      U(1,NY+1)=0.0D0+UB(1,NY+1)
      V(1,NY+1)=0.0D0+VB(1,NY+1)
      U(NX+1,NY+1)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,NY+1)-4.0D0*H(NX,NY+1)+
     13.0D0*H(NX+1,NY+1))*DXNEW)**2+((H(NX+1,1)-H(NX+1,NY))*DYNEW)**2)**
     2((XN-1.0D0)*0.5D0)*(H(NX+1,NY+1)**(XN+1.0D0))*XN2*DABS(H(NX-1,NY+1
     3)-4.0D0*H(NX,NY+1)+3.0D0*H(NX+1,NY+1))*DXNEW
     3+2.0D0*A*(DXX(NX+1,NY+1)**2+DXY(NX+1,NY+1)**2)*(-RHO*G)*H(NX+1,NY+
     41)**2*(H(NX-1,NY+1)-4.0D0*H(NX,NY+1)+3.0D0*H(NX+1,NY+1))*DXNEW/XN
     5+UB(NX+1,NY+1)
      V(NX+1,NY+1)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,NY+1)-4.0D0*H(NX,NY+1)+
     13.0D0*H(NX+1,NY+1))*DXNEW)**2+((H(NX+1,1)-H(NX+1,NY))*DYNEW)**2)**
     2((XN-1.0D0)*0.5D0)*H(NX+1,NY+1)**(XN+1.0D0)*XN2*(H(NX+1,1)-H(NX+1,
     3NY))*DYNEW
     3+2.0D0*A*(DXX(NX+1,NY+1)**2+DXY(NX+1,NY+1)**2)*(-RHO*G)*H(NX+1,NY+
     41)**2*(H(NX+1,1)-H(NX+1,NY))*DYNEW/XN+VB(NX+1,NY+1)
      RETURN
      END
      SUBROUTINE DIFFCOE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
C     DIMENSIONED COMMON BLOCKS
      COMMON /A/ H(200,200),U(200,200),V(200,200),D(200,200)
      COMMON /B/ OLDH(200,200),UNEW(200,200),VNEW(200,200)
      COMMON /C/ QNEWX(200,200),QNEWY(200,200),DELTHK(200,200),
     1           BOTTOM(200)
      COMMON /C3/ DXX(200,200),DXY(200,200)
      COMMON /C4/ OLDU(200,200),OLDV(200,200),OLDDXX(200,200),OLDDXY(200
     1,200),UB(200,200),VB(200,200)
C     PARAMETER COMMON BLOCKS
      COMMON /D/ BDOT,A,TB,XN,DT,RHO,G,DXL,DYL,SECYR
C     CALCULATED COMMON BLOCKS
      COMMON /E/ NX,NY,DXNEW,DYNEW,XN2,RG,X,DX,DY
C******MAIN BODY LOOPS*************************************************C
      DO 200 I=2,NX
      DO 201 J=2,NY
      D(I,J)=(2.0D0*A*(RHO*G)**(XN)*(H(I,J)**(XN+2.0D0))*
     1(((H(I+1,J)-H(I-1,J))*DXNEW)**2+((H(I,J+1)-H(I,J-1))*DYNEW
     2)**2)**((XN-1.0D0)*0.5D0))*XN2
     3+2.0D0*A*(DXX(I,J)**2+DXY(I,J)**2)*(-RHO*G)*H(I,J)**3
     4/XN
C     5+UB(I,J)/((H(I+1,J)-H(I-1,J))*DXNEW)
C     6+VB(I,J)/((H(I,J+1)-H(I,J-1))*DYNEW)
201   CONTINUE
200   CONTINUE
C*********ICE DIVIDE, END OF GLACIER LOOPS***************************C
      DO 202 J=2,NY
      D(1,J)=2.0D0*A*(RHO*G)**XN*H(1,J)**(XN+2.0D0)*(((H(1,J+1)-H(1,J-1
     1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*XN2
     3+2.0D0*A*(DXX(1,J)**2+DXY(1,J)**2)*(-RHO*G)*H(1,J)**3
     4/XN
C     6+VB(1,J)/((H(1,J+1)-H(1,J-1))*DYNEW)
      D(NX+1,J)=2.0D0*A*(RHO*G)**XN*H(NX+1,J)**(XN+2.0D0)*(((H(NX-1,J)-4
     1.0D0*H(NX,J)+3.0D0*H(NX+1,J))*DXNEW)**2+((H(NX+1,J+1)-H(NX+1,J-1))
     2*DYNEW)**2)**((XN-1.0D0)*0.5D0)*XN2
     3+2.0D0*A*(DXX(NX+1,J)**2+DXY(NX+1,J)**2)*(-RHO*G)*H(NX+1,J)**3
     4/XN
C     5+UB(NX+1,J)/((H(NX-1,J)-4.0D0*H(NX,J)+3.0D0*H(NX+1,J))*DXNEW)
C     6+VB(NX+1,J)/((H(NX+1,J+1)-H(NX+1,J-1))*DYNEW)
202   CONTINUE
C**********TOP,BOTTOM LOOPS******************************************C
      DO 203 I=2,NX
      D(I,1)=(2.0D0*A*(RHO*G)**(XN)*(H(I,1)**(XN+2.0D0))*
     1(((H(I+1,1)-H(I-1,1))*DXNEW)**2+((H(I,2)-H(I,NY+1))*DYNEW
     2)**2)**((XN-1.0D0)*0.5D0))*XN2
     3+2.0D0*A*(DXX(I,1)**2+DXY(I,1)**2)*(-RHO*G)*H(I,1)**3
     4/XN
C     5+UB(I,1)/((H(I+1,1)-H(I-1,1))*DXNEW)
C     6+VB(I,1)/((H(I,2)-H(I,NY+1))*DYNEW)
      D(I,NY+1)=(2.0D0*A*(RHO*G)**(XN)*(H(I,NY+1)**(XN+2.0D0))*
     1(((H(I+1,NY+1)-H(I-1,NY+1))*DXNEW)**2+((H(I,1)-H(I,NY))*DYNEW)
     2**2)**((XN-1.0D0)*0.5D0))*XN2
     3+2.0D0*A*(DXX(I,NY+1)**2+DXY(I,NY+1)**2)*(-RHO*G)*H(I,NY+1)**3
     4/XN
C     5+UB(I,NY+1)/((H(I+1,NY+1)-H(I-1,NY+1))*DXNEW)
C     6+VB(I,NY+1)/((H(I,1)-H(I,NY))*DYNEW)
203   CONTINUE
C**********THE 4 CORNERS*********************************************C
      D(1,1)=2.0D0*A*(RHO*G)**XN*H(1,1)**(XN+2.0D0)*(((H(1,2)-H(1,NY+1))
     1*DYNEW)**2)**((XN-1.0D0)*0.5D0)*XN2
     3+2.0D0*A*(DXX(1,1)**2+DXY(1,1)**2)*(-RHO*G)*H(1,1)**3
     4/XN
C     6+VB(1,1)/((H(1,2)-H(1,NY+1))*DYNEW)
      D(NX+1,1)=2.0D0*A*(RHO*G)**XN*H(NX+1,1)**(XN+2.0D0)*(((H(NX-1,1)-4
     1.0D0*H(NX,1)+3.0D0*H(NX+1,1))*DXNEW)**2+((H(NX+1,2)-H(NX+1,NY+1))*
     2DYNEW)**2)**((XN-1.0D0)*0.5D0)*XN2
     3+2.0D0*A*(DXX(NX+1,1)**2+DXY(NX+1,1)**2)*(-RHO*G)*H(NX+1,1)**3
     4/XN
C     5+UB(NX+1,1)/((H(NX-1,1)-4.0D0*H(NX,J)+3.0D0*H(NX+1,1))*DXNEW)
C     6+VB(NX+1,1)/((H(NX+1,2)-H(NX+1,NY+1))*DYNEW)
      D(1,NY+1)=2.0D0*A*(RHO*G)**XN*H(1,NY+1)**(XN+2.0D0)*(((H(1,1)-H(1,
     1NY))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*XN2
     3+2.0D0*A*(DXX(1,NY+1)**2+DXY(1,NY+1)**2)*(-RHO*G)*H(1,NY+1)**3
     4/XN
C     6+VB(1,NY+1)/((H(1,1)-H(1,NY))*DYNEW)
      D(NX+1,NY+1)=2.0D0*A*(RHO*G)**XN*H(NX+1,NY+1)**(XN+2.0D0)*(((H(NX-
     11,NY+1)-4.0D0*H(NX,NY+1)+3.0D0*H(NX+1,NY+1))*DXNEW)**2+((H(NX+1,1)
     2-H(NX+1,NY))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*XN2
     3+2.0D0*A*(DXX(NX+1,NY+1)**2+DXY(NX+1,NY+1)**2)*(-RHO*G)*H(NX+1,NY+
     41)**2/XN
C    5+UB(NX+1,NY+1)/((H(NX-1,NY+1)-4.0D0*H(NX,NY+1)+3.0D0*H(NX+1,NY+1))
C     6*DXNEW)+VB(NX+1,NY+1)/((H(NX+1,1)-H(NX+1,NY))*DYNEW)
C     6*DXNEW)
      RETURN
      END
      SUBROUTINE INIT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
C     DIMENSIONED COMMON BLOCKS
      COMMON /A/ H(200,200),U(200,200),V(200,200),D(200,200)
      COMMON /B/ OLDH(200,200),UNEW(200,200),VNEW(200,200)
      COMMON /C/ QNEWX(200,200),QNEWY(200,200),DELTHK(200,200),
     1           BOTTOM(200)
      COMMON /C3/ DXX(200,200),DXY(200,200)
      COMMON /C4/ OLDU(200,200),OLDV(200,200),OLDDXX(200,200),OLDDXY(200
     1,200),UB(200,200),VB(200,200)
C     PARAMETER COMMON BLOCKS
      COMMON /D/ BDOT,A,TB,XN,DT,RHO,G,DXL,DYL,SECYR
C     CALCULATED COMMON BLOCKS
      COMMON /E/ NX,NY,DXNEW,DYNEW,XN2,RG,X,DX,DY
      XXX=2.0D0*DX/10.0D0
      DO 101 J=1,NY+1
      X=0.0D0
      DO 10 I=1,NX+1
C      H(I,J)=DSQRT(2.0D0*TB*RG*(DXL+XXX-X))
C      H(I,J)=500.0D0
      H(I,J)=(1.0D0-(X/(DXL+XXX))**(4.D0/3.D0))**(3.D0/8.D0)*1963.32518
      OLDH(I,J)=H(I,J)
      X=X+DX
10    CONTINUE
101   CONTINUE
      DO 91 J=1,NY+1
      BOTTOM(J)=H(NX+1,J)
91    CONTINUE
      RETURN
      END
      SUBROUTINE HIEGHT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
C     DIMENSIONED COMMON BLOCKS
      COMMON /A/ H(200,200),U(200,200),V(200,200),D(200,200)
      COMMON /B/ OLDH(200,200),UNEW(200,200),VNEW(200,200)
      COMMON /C/ QNEWX(200,200),QNEWY(200,200),DELTHK(200,200),
     1           BOTTOM(200)
      COMMON /C3/ DXX(200,200),DXY(200,200)
      COMMON /C4/ OLDU(200,200),OLDV(200,200),OLDDXX(200,200),OLDDXY(200
     1,200),UB(200,200),VB(200,200)
C     PARAMETER COMMON BLOCKS
      COMMON /D/ BDOT,A,TB,XN,DT,RHO,G,DXL,DYL,SECYR
C     CALCULATED COMMON BLOCKS
      COMMON /E/ NX,NY,DXNEW,DYNEW,XN2,RG,X,DX,DY
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
      SUBROUTINE STRESS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CHGDXX(200,200),CHGDXY(200,200)
      COMMON /A/ H(200,200),U(200,200),V(200,200),D(200,200)
      COMMON /B/ OLDH(200,200),UNEW(200,200),VNEW(200,200)
      COMMON /C/ QNEWX(200,200),QNEWY(200,200),DELTHK(200,200),
     1           BOTTOM(200)
      COMMON /C3/ DXX(200,200),DXY(200,200)
      COMMON /C4/ OLDU(200,200),OLDV(200,200),OLDDXX(200,200),OLDDXY(200
     1,200),UB(200,200),VB(200,200)
C     PARAMETER COMMON BLOCKS
      COMMON /D/ BDOT,A,TB,XN,DT,RHO,G,DXL,DYL,SECYR
C     CALCULATED COMMON BLOCKS
      COMMON /E/ NX,NY,DXNEW,DYNEW,XN2,RG,X,DX,DY
      DO 21 I=1,NX+1
      DO 20 J=1,NY+1
      OLDDXX(I,J)=DXX(I,J)
      OLDDXY(I,J)=DXY(I,J)
20    CONTINUE
21    CONTINUE
C******MAIN BODY LOOP**************************************************C
      DO 101 I=2,NX
      DO 100 J=2,NY
      DXX(I,J)=(U(I+1,J)-U(I-1,J))*DXNEW/(A*(DXX(I,J)**2+DXY(I,J)**2
     1+(RHO*G*H(I,J))**2*(((H(I+1,J)-H(I-1,J))*DXNEW)**2+
     3((H(I,J+1)-H(I,J-1))*DYNEW)**2)))
      DXY(I,J)=0.5D0*((U(I,J+1)-U(I,J-1))*DYNEW+(V(I+1,J)-V(I-1,J))
     1*DXNEW)/(A*(DXX(I,J)**2+DXY(I,J)**2
     2+(RHO*G*H(I,J))**2*(((H(I+1,J)-H(I-1,J))*DXNEW)**2+
     3((H(I,J+1)-H(I,J-1))*DYNEW)**2)))
100   CONTINUE
101   CONTINUE
C******ICE DIVIDE, END OF GLACIER LOOP*********************************C
      DO 102 J=2,NY
      DXX(1,J)=2.0D0*U(2,J)*DXNEW/(A*(DXX(1,J)**2+DXY(1,J)**2
     1+(RHO*G*H(1,J))**2*(
     3((H(1,J+1)-H(1,J-1))*DYNEW)**2)))
      DXY(1,J)=0.5D0*((U(1,J+1)-U(1,J-1))*DYNEW+(2.0D0*V(2,J))
     1*DXNEW)/(A*(DXX(1,J)**2+DXY(1,J)**2
     2+(RHO*G*H(1,J))**2*(
     3((H(I,J+1)-H(I,J-1))*DYNEW)**2)))
      DXX(NX+1,J)=(U(NX-1,J)-4.0D0*U(NX,J)+3.0D0*U(NX+1,J))*DXNEW/(A*(DX
     1X(NX+1,J)**2+DXY(NX+1,J)**2+(RHO*G*H(NX+1,J))**2*(((H(NX-1,J)-4.0D
     30*H(NX,J)+3.0D0*H(NX+1,J))*DXNEW)**2+((H(NX+1,J+1)-H(NX+1,J-1))*DY
     4NEW)**2)))
      DXY(NX+1,J)=0.5D0*((U(NX+1,J+1)-U(NX+1,J-1))*DYNEW+(V(NX-1,J)-4.0D
     10*V(NX,J)+3.0D0*V(NX+1,J))*DXNEW)/(A*(DXX(NX+1,J)**2+DXY(NX+1,J)**
     22+(RHO*G*H(NX+1,J))**2*(((H(NX-1,J)-4.0D0*H(NX,J)+3.0D0*H(NX+1,J))
     3*DXNEW)**2+((H(NX+1,J+1)-H(NX+1,J-1))*DYNEW)**2)))
102   CONTINUE
C******TOP AND BOTTOM LOOPS (THE JOIN)*********************************C
      DO 103 I=2,NX
      DXX(I,1)=(U(I+1,1)-U(I-1,1))*DXNEW/(A*(DXX(I,1)**2+DXY(I,1)**2
     1+(RHO*G*H(I,1))**2*(((H(I+1,1)-H(I-1,1))*DXNEW)**2+
     3((H(I,2)-H(I,NY+1))*DYNEW)**2)))
      DXY(I,1)=0.5D0*((U(I,2)-U(I,NY+1))*DYNEW+(V(I+1,1)-V(I-1,1))
     1*DXNEW)/(A*(DXX(I,1)**2+DXY(I,1)**2
     2+(RHO*G*H(I,1))**2*(((H(I+1,1)-H(I-1,1))*DXNEW)**2+
     3((H(I,2)-H(I,NY+1))*DYNEW)**2)))
      DXX(I,NY+1)=(U(I+1,NY+1)-U(I-1,NY+1))*DXNEW/(A*(DXX(I,NY+1)**2+DXY
     1(I,NY+1)**2+(RHO*G*H(I,NY+1))**2*(((H(I+1,NY+1)-H(I-1,NY+1))*DXNEW
     3)**2+((H(I,2)-H(I,NY))*DYNEW)**2)))
      DXY(I,NY+1)=0.5D0*((U(I,1)-U(I,NY))*DYNEW+(V(I+1,NY+1)-V(I-1,NY+1)
     1)*DXNEW)/(A*(DXX(I,NY+1)**2+DXY(I,NY+1)**2
     2+(RHO*G*H(I,NY+1))**2*(((H(I+1,NY+1)-H(I-1,NY+1))*DXNEW)**2+
     3((H(I,1)-H(I,NY))*DYNEW)**2)))
103   CONTINUE
C******THE 4 CORNERS***************************************************C
      DXX(1,1)=(2.0D0*U(2,J))*DXNEW/(A*(DXX(1,1)**2+DXY(1,1)**2
     1+(RHO*G*H(1,1))**2*(
     3((H(1,2)-H(1,NY+1))*DYNEW)**2)))
      DXY(1,1)=0.5D0*((U(1,2)-U(1,NY+1))*DYNEW+(2.0D0*V(2,1))
     1*DXNEW)/(A*(DXX(1,1)**2+DXY(1,1)**2
     2+(RHO*G*H(1,1))**2*(
     3((H(1,2)-H(1,NY+1))*DYNEW)**2)))
      DXX(1,NY+1)=(2.0D0*U(2,NY+1))*DXNEW/(A*(DXX(1,NY+1)**2+DXY(1,NY+1)
     1**2+(RHO*G*H(1,NY+1))**2*(
     3((H(1,1)-H(1,NY))*DYNEW)**2)))
      DXY(1,NY+1)=0.5D0*((U(1,1)-U(1,NY))*DYNEW+(2.0D0*V(2,NY+1))
     1*DXNEW)/(A*(DXX(1,NY+1)**2+DXY(1,NY+1)**2
     2+(RHO*G*H(1,NY+1))**2*(
     3((H(1,1)-H(1,NY))*DYNEW)**2)))
      DXX(NX+1,1)=(U(NX-1,1)-4.0D0*U(NX,1)+3.0D0*U(NX+1,1))*DXNEW/(A*(DX
     1X(NX+1,1)**2+DXY(NX+1,1)**2+(RHO*G*H(NX+1,1))**2*(((H(NX-1,1)-4.0D
     30*H(NX,1)+3.0D0*H(NX+1,1))*DXNEW)**2+((H(NX+1,2)-H(NX+1,NY+1))*DYN
     4EW)**2)))
      DXY(NX+1,1)=0.5D0*((U(NX+1,2)-U(NX+1,NY+1))*DYNEW+(V(NX-1,1)-4.0D0
     1*V(NX,1)+3.0D0*V(NX+1,1))*DXNEW)/(A*(DXX(NX+1,1)**2+DXY(NX+1,1)**2
     2+(RHO*G*H(NX+1,1))**2*(((H(NX-1,1)-4.0D0*H(NX,1)+3.0D0*H(NX+1,1))*
     3DXNEW)**2+((H(NX+1,2)-H(NX+1,NY+1))*DYNEW)**2)))
      DXX(NX+1,NY+1)=(U(NX-1,NY+1)-4.0D0*U(NX,NY+1)+3.0D0*U(NX+1,NY+1))*
     1DXNEW/(A*(DXX(NX+1,NY+1)**2+DXY(NX+1,NY+1)**2+(RHO*G*H(NX+1,NY+1))
     3**2*(((H(NX-1,NY+1)-4.0D0*H(NX,NY+1)+3.0D0*H(NX+1,NY+1))*DXNEW)**2
     4+((H(NX+1,1)-H(NX+1,NY))*DYNEW)**2)))
      DXY(NX+1,NY+1)=0.5D0*((U(NX+1,1)-U(NX+1,NY))*DYNEW+(V(NX-1,NY+1)-4
     1.0D0*V(NX,NY+1)+3.0D0*V(NX+1,NY+1))*DXNEW)/(A*(DXX(NX+1,NY+1)**2+D
     2XY(NX+1,NY+1)**2+(RHO*G*H(NX+1,NY+1))**2*(((H(NX-1,NY+1)-4.0D0*H(N
     3X,NY+1)+3.0D0*H(NX+1,NY+1))*DXNEW)**2+((H(NX+1,1)-H(NX+1,NY))*DYNE
     4W)**2)))
      RETURN
      END
      SUBROUTINE CHK(NX,NY)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      WRITE(6,200) NX,NY
200   FORMAT(' ','NX= ',I4,5X,'NY= ',I4)
      RETURN
      END
      SUBROUTINE VELINT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CHGU(200,200),CHGV(200,200)
C     DIMENSIONED COMMON BLOCKS
      COMMON /A/ H(200,200),U(200,200),V(200,200),D(200,200)
      COMMON /B/ OLDH(200,200),UNEW(200,200),VNEW(200,200)
      COMMON /C/ QNEWX(200,200),QNEWY(200,200),DELTHK(200,200),
     1           BOTTOM(200)
      COMMON /C3/ DXX(200,200),DXY(200,200)
      COMMON /C4/ OLDU(200,200),OLDV(200,200),OLDDXX(200,200),OLDDXY(200
     1,200),UB(200,200),VB(200,200)
C     PARAMETER COMMON BLOCKS
      COMMON /D/ BDOT,A,TB,XN,DT,RHO,G,DXL,DYL,SECYR
C     CALCULATED COMMON BLOCKS
      COMMON /E/ NX,NY,DXNEW,DYNEW,XN2,RG,X,DX,DY
      DO 300 I=1,NX+1
      DO 301 J=1,NY+1
      OLDU(I,J)=U(I,J)
      OLDV(I,J)=V(I,J)
301   CONTINUE
300   CONTINUE
C******MAIN BODY LOOPS*************************************************C
      DO 200 I=2,NX
      DO 201 J=2,NY
      U(I,J)=2.0D0*A*(RHO*G)**XN*(((H(I+1,J)-H(I-1,J))*DXNEW)**2+((
     1H(I,J+1)-H(I,J-1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,J)**(
     2XN+1.0D0))*XN2*DABS(H(I+1,J)-H(I-1,J))*DXNEW
     3+2.0D0*A*(DXX(I,J)**2+DXY(I,J)**2)*(-RHO*G)*H(I,J)**2*(H(I+1,J)-H(
     4I-1,J))*DXNEW/XN
      V(I,J)=2.0D0*A*(RHO*G)**XN*(((H(I+1,J)-H(I-1,J))*DXNEW)**2+((
     1H(I,J+1)-H(I,J-1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,J)**(
     2XN+1.0D0))*XN2*(H(I,J+1)-H(I,J-1))*DYNEW
     3+2.0D0*A*(DXX(I,J)**2+DXY(I,J)**2)*(-RHO*G)*H(I,J)**2*(H(I,J+1)-H(
     4I,J-1))*DYNEW/XN
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
     3+2.0D0*A*(DXX(NX+1,J)**2+DXY(NX+1,J)**2)*(-RHO*G)*H(NX+1,J)**2*(H(
     4NX-1,J)-4.0D0*H(NX,J)+3.0D0*H(NX+1,J))*DXNEW/XN
      V(NX+1,J)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,J)-4.0D0*H(NX,J)+3.0D0*H(N
     1X+1,J))*DXNEW)**2+((H(NX+1,J+1)-H(NX+1,NY+1))*DYNEW)**2)**((XN-1.0
     2D0)*0.5D0)*H(NX+1,J)**(XN+1.0D0)*XN2*(H(NX+1,J+1)-H(NX+1,J-1))*DYN
     3EW
     3+2.0D0*A*(DXX(NX+1,J)**2+DXY(NX+1,J)**2)*(-RHO*G)*H(NX+1,J)**2*(H(
     4NX+1,J+1)-H(NX+1,J-1))*DYNEW/XN
202   CONTINUE
C**********TOP,BOTTOM LOOPS******************************************C
      DO 203 I=2,NX
      U(I,1)=2.0D0*A*(RHO*G)**XN*(((H(I+1,1)-H(I-1,1))*DXNEW)**2+((
     1H(I,2)-H(I,NY+1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,1)**(
     2XN+1.0D0))*XN2*DABS(H(I+1,1)-H(I-1,1))*DXNEW
     3+2.0D0*A*(DXX(I,1)**2+DXY(I,1)**2)*(-RHO*G)*H(I,1)**2*(H(I+1,1)-H(
     4I-1,1))*DXNEW/XN
      V(I,1)=2.0D0*A*(RHO*G)**XN*(((H(I+1,1)-H(I-1,1))*DXNEW)**2+((
     1H(I,2)-H(I,NY+1))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,1)**(XN+1
     2.0D0))*XN2*(H(I,2)-H(I,NY+1))*DYNEW
     3+2.0D0*A*(DXX(I,1)**2+DXY(I,1)**2)*(-RHO*G)*H(I,1)**2*(H(I,2)-H(
     4I,NY+1))*DYNEW/XN
      U(I,NY+1)=2.0D0*A*(RHO*G)**XN*(((H(I+1,NY+1)-H(I-1,NY+1))*DXNEW)**
     12+((H(I,1)-H(I,NY))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,NY+1)**(XN
     2+1.0D0))*XN2*DABS(H(I+1,NY+1)-H(I-1,NY+1))*DXNEW
     3+2.0D0*A*(DXX(I,NY+1)**2+DXY(I,NY+1)**2)*(-RHO*G)*H(I,NY+1)**2*(H(
     4I+1,NY+1)-H(I-1,NY+1))*DXNEW/XN
      V(I,NY+1)=2.0D0*A*(RHO*G)**XN*(((H(I+1,NY+1)-H(I-1,NY+1))*DXNEW)**
     12+((H(I,1)-H(I,NY))*DYNEW)**2)**((XN-1.0D0)*0.5D0)*(H(I,NY+1)**(XN
     2+1.0D0))*XN2*(H(I,1)-H(I,NY))*DYNEW
     3+2.0D0*A*(DXX(I,NY+1)**2+DXY(I,NY+1)**2)*(-RHO*G)*H(I,NY+1)**2*(H(
     4I,1)-H(I,NY))*DYNEW/XN
203   CONTINUE
C**********THE 4 CORNERS*********************************************C
      U(1,1)=0.0D0
      V(1,1)=0.0D0
      U(NX+1,1)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,1)-4.0D0*H(NX,1)+3.0D0*H(N
     1X+1,1))*DXNEW)**2+((H(NX+1,2)-H(NX+1,NY+1))*DYNEW)**2)**((XN-1.0D
     20)*0.5D0)*(H(NX+1,1)**(XN+1.0D0))*XN2*DABS(H(NX-1,1)-4.0D0*H(NX,
     31)+3.0D0*H(NX+1,1))*DXNEW
     3+2.0D0*A*(DXX(NX+1,1)**2+DXY(NX+1,1)**2)*(-RHO*G)*H(NX+1,1)**2*(H(
     4NX-1,1)-4.0D0*H(NX,1)+3.0D0*H(NX+1,1))*DXNEW/XN
      V(NX+1,1)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,1)-4.0D0*H(NX,1)+3.0D0*H(N
     1X+1,1))*DXNEW)**2+((H(NX+1,2)-H(NX+1,NY+1))*DYNEW)**2)**((XN-1.0D
     20)*0.5D0)*H(NX+1,1)**(XN+1.0D0)*XN2*(H(NX+1,2)-H(NX+1,NY+1))*DYNEW
     3+2.0D0*A*(DXX(NX+1,1)**2+DXY(NX+1,1)**2)*(-RHO*G)*H(NX+1,1)**2*(H(
     4NX+1,2)-H(NX+1,NY+1))*DYNEW/XN
      U(1,NY+1)=0.0D0
      V(1,NY+1)=0.0D0
      U(NX+1,NY+1)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,NY+1)-4.0D0*H(NX,NY+1)+
     13.0D0*H(NX+1,NY+1))*DXNEW)**2+((H(NX+1,1)-H(NX+1,NY))*DYNEW)**2)**
     2((XN-1.0D0)*0.5D0)*(H(NX+1,NY+1)**(XN+1.0D0))*XN2*DABS(H(NX-1,NY+1
     3)-4.0D0*H(NX,NY+1)+3.0D0*H(NX+1,NY+1))*DXNEW
     3+2.0D0*A*(DXX(NX+1,NY+1)**2+DXY(NX+1,NY+1)**2)*(-RHO*G)*H(NX+1,NY+
     41)**2*(H(NX-1,NY+1)-4.0D0*H(NX,NY+1)+3.0D0*H(NX+1,NY+1))*DXNEW/XN
      V(NX+1,NY+1)=2.0D0*A*(RHO*G)**XN*(((H(NX-1,NY+1)-4.0D0*H(NX,NY+1)+
     13.0D0*H(NX+1,NY+1))*DXNEW)**2+((H(NX+1,1)-H(NX+1,NY))*DYNEW)**2)**
     2((XN-1.0D0)*0.5D0)*H(NX+1,NY+1)**(XN+1.0D0)*XN2*(H(NX+1,1)-H(NX+1,
     3NY))*DYNEW
     3+2.0D0*A*(DXX(NX+1,NY+1)**2+DXY(NX+1,NY+1)**2)*(-RHO*G)*H(NX+1,NY+
     41)**2*(H(NX+1,1)-H(NX+1,NY))*DYNEW/XN
      DO 21 I=1,NX+1
      DO 25 J=1,NY+1
      UNEW(I,J)=U(I,J)*SECYR
      VNEW(I,J)=V(I,J)*SECYR
      QNEWX(I,J)=U(I,J)*H(I,J)*(1.0D0)*SECYR
      QNEWY(I,J)=V(I,J)*H(I,J)*(1.0D0)*SECYR
25    CONTINUE
21    CONTINUE
      RETURN
      END
      SUBROUTINE SIGINT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CHGDXX(200,200),CHGDXY(200,200)
      COMMON /A/ H(200,200),U(200,200),V(200,200),D(200,200)
      COMMON /B/ OLDH(200,200),UNEW(200,200),VNEW(200,200)
      COMMON /C/ QNEWX(200,200),QNEWY(200,200),DELTHK(200,200),
     1           BOTTOM(200)
      COMMON /C3/ DXX(200,200),DXY(200,200)
      COMMON /C4/ OLDU(200,200),OLDV(200,200),OLDDXX(200,200),OLDDXY(200
     1,200),UB(200,200),VB(200,200)
C     PARAMETER COMMON BLOCKS
      COMMON /D/ BDOT,A,TB,XN,DT,RHO,G,DXL,DYL,SECYR
C     CALCULATED COMMON BLOCKS
      COMMON /E/ NX,NY,DXNEW,DYNEW,XN2,RG,X,DX,DY
      DO 21 I=1,NX+1
      DO 20 J=1,NY+1
      OLDDXX(I,J)=DXX(I,J)
      OLDDXY(I,J)=DXY(I,J)
20    CONTINUE
21    CONTINUE
C******MAIN BODY LOOP**************************************************C
      DO 101 I=2,NX
      DO 100 J=2,NY
      DXX(I,J)=(DABS(U(I+1,J)-U(I-1,J))*DXNEW/A)**(1.D0/3.D0)
      DXY(I,J)=(0.5D0*(DABS(U(I,J+1)-U(I,J-1))*DYNEW+DABS(V(I+1,J)-V(I-1
     1,J))*DXNEW)/A)**(1.D0/3.D0)
100   CONTINUE
101   CONTINUE
C******ICE DIVIDE, END OF GLACIER LOOP*********************************C
      DO 102 J=2,NY
      DXX(1,J)=(2.0D0*U(2,J)*DXNEW/A)**(1.D0/3.D0)
      DXY(1,J)=(0.5D0*(DABS(U(1,J+1)-U(1,J-1))*DYNEW+(2.0D0*V(2,J))
     1*DXNEW)/A)**(1.D0/3.D0)
      DXX(NX+1,J)=(DABS(U(NX-1,J)-4.0D0*U(NX,J)+3.0D0*U(NX+1,J))*DXNEW/A
     1)**(1.D0/3.D0)
      DXY(NX+1,J)=(0.5D0*(DABS(U(NX+1,J+1)-U(NX+1,J-1))*DYNEW+DABS(V(NX-
     11,J)-4.0D0*V(NX,J)+3.0D0*V(NX+1,J))*DXNEW)/A)**(1.D0/3.D0)
102   CONTINUE
C******TOP AND BOTTOM LOOPS (THE JOIN)*********************************C
      DO 103 I=2,NX
      DXX(I,1)=(DABS(U(I+1,1)-U(I-1,1))*DXNEW/A)**(1.D0/3.D0)
      DXY(I,1)=(0.5D0*(DABS(U(I,2)-U(I,NY+1))*DYNEW+DABS(V(I+1,1)-V(I-1,
     11))*DXNEW)/A)**(1.D0/3.D0)
      DXX(I,NY+1)=(DABS(U(I+1,NY+1)-U(I-1,NY+1))*DXNEW/A)**(1.D0/3.D0)
      DXY(I,NY+1)=(0.5D0*(DABS(U(I,1)-U(I,NY))*DYNEW+DABS(V(I+1,NY+1)-V(
     1I-1,NY+1))*DXNEW)/A)**(1.D0/3.D0)
103   CONTINUE
C******THE 4 CORNERS***************************************************C
      DXX(1,1)=((2.0D0*U(2,J))*DXNEW/A)**(1.D0/3.D0)
      DXY(1,1)=(0.5D0*(DABS(U(1,2)-U(1,NY+1))*DYNEW+(2.0D0*V(2,1))
     1*DXNEW)/A)**(1.D0/3.D0)
      DXX(1,NY+1)=((2.0D0*U(2,NY+1))*DXNEW/A)**(1.D0/3.D0)
      DXY(1,NY+1)=(0.5D0*(DABS(U(1,1)-U(1,NY))*DYNEW+(2.0D0*V(2,NY+1))
     1*DXNEW)/A)**(1.D0/3.D0)
      DXX(NX+1,1)=(DABS(U(NX-1,1)-4.0D0*U(NX,1)+3.0D0*U(NX+1,1))*DXNEW/A
     1)**(1.D0/3.D0)
      DXY(NX+1,1)=(0.5D0*(DABS(U(NX+1,2)-U(NX+1,NY+1))*DYNEW+DABS(V(NX-1
     1,1)-4.0D0*V(NX,1)+3.0D0*V(NX+1,1))*DXNEW)/A)**(1.D0/3.D0)
      DXX(NX+1,NY+1)=(DABS(U(NX-1,NY+1)-4.0D0*U(NX,NY+1)+3.0D0*U(NX+1,NY
     1+1))*DXNEW/A)**(1.D0/3.D0)
      DXY(NX+1,NY+1)=(0.5D0*(DABS(U(NX+1,1)-U(NX+1,NY))*DYNEW+DABS(V(NX-
     11,NY+1)-4.0D0*V(NX,NY+1)+3.0D0*V(NX+1,NY+1))*DXNEW)/A)**(1.D0/3.D0
     2)
      RETURN
      END
//GO.SYSIN DD *
/*
//
