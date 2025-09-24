C---------------------------------------------------------------------------------
      PROGRAM const_reg_gene_exp

c     integrates NEQU first order explicit ODE's using RK4 algorithm
c     The equations to be integrated are to be Defined in subroutine DERIVADES
c     RK4 is implemented in a subroutine and is independent of the specific ODEs
  
      IMPLICIT NONE
      INTEGER NEQU,activ,nwrite,write,i
      PARAMETER (NEQU=18,activ=4)
c     t: independent variable
c     h: integration step
c     y: vector of NEQU variables at independent variable value t
c     yout: vector of NEQU variables at independent variable value t+h
c     tf: final value of t to integrate (initial value is assumed to be 0)
    
      DOUBLE PRECISION y(nequ),yout(nequ),e(activ)
      DOUBLE PRECISION t,h,tf,e0(activ),step

c     parameters 
      DOUBLE PRECISION Kf, Kb, Kcat,e01Max

      COMMON /param/ Kf,Kb,Kcat,tf,h,write,nwrite
      
      OPEN(unit=20,file="Act_Kinases_concentrationsVSTime.dat",status="unknown")
      OPEN(unit=21,file="Fig1_stationaryConc_vs_EC50.dat",status="unknown")
      OPEN(unit=22,file="Stationary_concentrationVSE1_log.dat",status="unknown")


ccccccccccccccccccc parameter values cccccccccccccccccc

cccc  Here the reaction constants ar set up   ccccccccc
 
      Kf = 1.d0/150.d0
      Kb = 1.d0
      Kcat = 1.d0

ccccccccccccccccccc initial conditions ccccccccccccccccc      
      t=0.0d0

C y1 is MAPKKK, y5 MAPKK and y12 MAPK
      DO i = 1, nequ
        y(i)=0.d0
      ENDDO
      y(1)=3.d0
      y(5)=1200.d0
      y(12)=1200.d0

C Activators Ei
      e0(1)=1.d0
      e0(2)=0.3d0
      e0(3)=0.3d0
      e0(4)=120.d0

C It is necessary to save the initial values of Ei
      e(1)=e0(1)
      e(2)=e0(2)
      e(3)=e0(3)
      e(4)=e0(4)
      
cccccccccccccccccc integration over time cccccccccccccccc               
c until final time tf with time step h, saving data every nwrite cccccc      
   
      tf=10000.d0
      h=0.1d0
      nwrite=200

      write=1 !This will act as a boolean to know if we write the data simulation
      CALL DeterSimulation(activ,nequ,y,yout,e,e0)

C Let's see the acivated kinases dependance on E1
      e0(1)=0.d0

      e(1)=e0(1)
      e(2)=e0(2)
      e(3)=e0(3)
      e(4)=e0(4)

      write=0 !No need to write this data
      step=0.0005d0 !Increasing of E1_0
      e01Max=5d0 !We will study the simulation from 0 to 5

      WRITE(21,*)"#      EC50_KKK-P                    [KKK-P]",
     & "              EC50_KK-PP                    [KK-PP]",
     & "               EC50_K-PP                    [K-PP]"

      WRITE(22,*)"#        E1/1000                    [KKK-P]",
     & "                    E1/1000                  [KK-PP]",
     & "                    E1/1000                   [K-PP]"
      DO WHILE ((10.d0**(-10.d0)).le.(e01Max-e0(1)))
C Again the initial conditions
      DO i = 1, nequ
        y(i)=0.d0
      ENDDO
      y(1)=3.d0
      y(5)=1200.d0
      y(12)=1200.d0
C---------------------------------------------------------------------------------
      WRITE(20,*)"#          time               [KKK-P]               [KK-PP]",
     & "              [K-PP]"
      CALL DeterSimulation(activ,nequ,y,yout,e,e0)

C We will write the normalized concentration variation of the activated kinases vs the E1 in multiples of EC50
      WRITE(21,*)e0(1)/0.2974,yout(3)/2.75442258,e0(1)/0.01884d0,
     & yout(10)/1081.2175031,e0(1)/0.00674,yout(17)/989.0932284
C This is the normalized concentration variation of the activated kinases vs the E1 for a log graphic
      WRITE(22,*)e0(1)/1000.,yout(3)/2.75442258,e0(1)/1000.,
     & yout(10)/1081.2175031,e0(1)/1000.,yout(17)/989.0932284

      e0(1)=e0(1)+step !Update of E1_0

      e(1)=e0(1)
      e(2)=e0(2)
      e(3)=e0(3)
      e(4)=e0(4)
      ENDDO
C Here satures the activated kinases concentration for big E1 and time
C Y3 2.7544225763424928, y10 1081.2175031239531, y17 989.09322839745698
      WRITE(*,*)e0(1),yout(3)/2.7544225763424928,
     & yout(10)/1081.2175031239531,yout(17)/989.09322839745698

      CLOSE(20)
      CLOSE(21)
      CLOSE(22)

      END 

C-----------------------------------------------------------------------
C-----------------------------SUBROUTINES-------------------------------
C-----------------------------------------------------------------------
      subroutine DFUNCTIONS(activ,NEQU,t,y,dy,e)
      implicit none
      double precision y(nequ),dy(nequ),e(activ)
      double precision t,tf,h
      integer nequ,activ,write,nwrite

c     parameters 
      double precision Kf, Kb, Kcat

      common /param/ Kf,Kb,Kcat,tf,h,write,nwrite

C Bunch of differential equations which drives the dynamics
        dy(1)=-Kf*y(1)*e(1)+Kb*y(2)+Kcat*y(4)
        dy(2)=Kf*y(1)*e(1)-(Kb+Kcat)*y(2)
        dy(3)=Kcat*(y(2)+y(9)+y(6))-Kf*(y(3)*e(2)+y(5)*y(3)+y(7)*y(3))+
     &    Kb*(y(4)+y(6)+y(9))
        dy(4)=Kf*y(3)*e(2)-(Kcat+Kb)*y(4)
        dy(5)=-Kf*y(5)*y(3)+Kcat*y(8)+Kb*y(6)
        dy(6)=Kf*y(5)*y(3)-(Kcat+Kb)*y(6)

        dy(7)=Kcat*(y(6)+y(11))-Kf*(y(7)*e(3)+y(7)*y(3))+
     &    Kb*(y(8)+y(9))
        dy(8)=Kf*y(7)*e(3)-(Kb+Kcat)*y(8)
        dy(9)=Kf*y(7)*y(3)-(Kb+Kcat)*y(9)
        dy(10)=Kcat*(y(9)+y(16)+y(13))+Kb*(y(11)+y(13)+y(16))-
     &    Kf*(y(10)*e(3)+y(10)*y(12)+y(10)*y(14))
        dy(11)=Kf*y(10)*e(3)-(Kb+Kcat)*y(11)
        dy(12)=-Kf*y(12)*y(10)+Kb*y(13)+Kcat*y(15)

        dy(13)=Kf*y(10)*y(12)-(Kcat+Kb)*y(13)
        dy(14)=Kcat*(y(13)+y(18))+Kb*(y(15)+y(16))-
     &    Kf*(y(14)*e(4)+y(14)*y(10))
        dy(15)=Kf*y(14)*e(4)-(Kb+Kcat)*y(15)
        dy(16)=Kf*y(14)*y(10)-(Kb+Kcat)*y(16)
        dy(17)=Kcat*y(16)-Kf*y(17)*e(4)+Kb*y(18)
        dy(18)=Kf*y(17)*e(4)-(Kb+Kcat)*y(18)
      return
      end

C-----------------------------------------------------------------------
      subroutine RK4(activ,nequ,x,h,y,yout,e)
      implicit none 
      double precision y(nequ),yout(nequ)
      double precision dy1(nequ),dy2(nequ),dy3(nequ),dy4(nequ)
      double precision x1,x3,e(activ)
      double precision y1(nequ),y2(nequ),y3(nequ)
      double precision x,h
      integer nequ,activ
C The rungekutta4 methor for the integration of the ODEs
         call DFUNCTIONS(activ,NEQU,x,y,dy1,e)
         x1=x+h*0.5d0
         y1=y+h*0.5d0*dy1  
         call DFUNCTIONS(activ,NEQU,x1,y1,dy2,e)
         y2=y+h*0.5d0*dy2  
         call DFUNCTIONS(activ,NEQU,x1,y2,dy3,e)
         x3=x+h
         y3=y+h*dy3
         call DFUNCTIONS(activ,NEQU,x3,y3,dy4,e)

         yout=y+h*(dy1+2.0d0*(dy2+dy3)+dy4)/6.0d0
      return
      end

C-------------------------------------------------------------
      SUBROUTINE DeterSimulation(activ,nequ,y,yout,e,e0)
      IMPLICIT NONE
      INTEGER j,npas,nequ,activ,nwrite,i,write
      DOUBLE PRECISION y(nequ),yout(nequ),e(activ),e0(activ),
     & t,h,tf,Kf,Kb,Kcat

      COMMON /param/ Kf,Kb,Kcat,tf,h,write,nwrite
C This is for a whole simulation from t=0 to a t=tf for the MAPK cascade

      npas=int(tf/h)

      do j=1,npas
           
         CALL RK4(activ,nequ,t,h,y,yout,e)
ccccc    update independent variable
         t=dble(j)*h

ccccc    write data 
         IF((write.eq.1).and.(mod(j,nwrite).eq.0))then
         WRITE(20,200)t,yout(3),yout(10),yout(17)
         ENDIF
     
ccccc    update variables for next step
         DO i=1,nequ
            y(i)=yout(i)
         ENDDO

C Mass conservation to know the stimuli concentration at real time
         e(1)=e0(1)-y(2)
         e(2)=e0(2)-y(4)
         e(3)=e0(3)-y(8)-y(11)
         e(4)=e0(4)-y(15)-y(18)

      ENDDO

      IF (write.eq.1) THEN
      WRITE(20,*)
      WRITE(20,*)
      ENDIF 
    
200   FORMAT(5(e20.8,1x))

      END
