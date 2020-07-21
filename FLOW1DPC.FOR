       PROGRAM FLOW1D
       DIMENSION ELE(50),X(50),H(50),DELH(50)
       RHO=910.0
       G=9.8
       H(1)=100.0
       DELX=50000.0
       DELH(1)=0.00
       HT=100.00
       WRITE(0,10)
10     FORMAT(' ','TYPE IN THE BASAL DRAG')
       READ(*,11) TAU
11     FORMAT(F10.3)
       WRITE(0,12) 
12     FORMAT(' ','TYPE IN THE NUMBER OF ELEVATIONS')
       READ(*,13) N
13     FORMAT(I2)
       WRITE(0,14)
14     FORMAT(' ','TYPE IN THE ELEVATIONS ABOVE SEA LEVEL')
       READ(*,15) (ELE(I),I=1,N)
15     FORMAT(F10.3)
       HT=HT+ELE(1)
       DO 16 I=1,N
       DELH(I+1)=TAU*DELX/(RHO*G*H(I))
       H(I+1)=H(I)-(ELE(I+1)-ELE(I))+DELH(I+1)
       WRITE(0,17) I,DELX,H(I),DELH(I),HT
17     FORMAT(' ','I= ',I2,1X,'DELX= ',F10.3,1X,'H= ',F10.3,1X,'DELH= 
     1',F10.3,1X,'HT= ',F10.3)
       HT=HT+DELH(I+1)
16     CONTINUE
       END
