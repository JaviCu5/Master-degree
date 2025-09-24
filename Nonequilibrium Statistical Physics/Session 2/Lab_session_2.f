C23456789012345678901234567890123456789012345678901234567890123456789012

c-----PRÀCTICA 2: Monte Carlo simulations of Hard Disks----------

c-----Javier Castillo Uviña


c-----MAIN PROGRAM------------------------------------------------------
      PROGRAM project2
c-----Variables---------------------------------------------------------
      IMPLICIT NONE
      !SECTION 2
      INTEGER ISeed, NPart, tau
      DOUBLE PRECISION Sigma, L, ComputesL, Phis(3), Phi
      PARAMETER (tau = 5000, NPart = 200, ISeed = 16866231)

      CALL SRand(ISeed)

      Sigma = 1.d0
      Phis = (/ 0.05d0, 0.2d0, 0.5d0 /)
c-----------------------------------------------------------------------
      !EXERCICE 1 
c-----------------------------------------------------------------------
      Phi = 0.3d0
      L = ComputesL(NPart, Phi, sigma)
      WRITE(*,*)L

C     CALL Excercise1(NPart, Sigma, L, tau)

c-----------------------------------------------------------------------
      !EXERCICE 2 
c-----------------------------------------------------------------------

C      CALL Excercise2(NPart, Sigma, tau)


c-----------------------------------------------------------------------
      !EXERCICE 3, 4 and 5
c-----------------------------------------------------------------------
      CALL Excercise5(NPart, Sigma, Phi, tau)


C      CALL system ("gnuplot -p Grafics2.gnu")
      CALL system ("gnuplot -p Grafics2Gravity.gnu")

      END PROGRAM project2
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                            END PROGRAM
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                             FUNCTIONS
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

C-----Periodic Boundary Condition (PBC)---------------------------------
      DOUBLE PRECISION FUNCTION PBCCheck(Var, Increment, L)
       IMPLICIT NONE 
       DOUBLE PRECISION Var, Increment, L, NewVar

       NewVar = Var + Increment
        IF (NewVar.GT.L) THEN
          PBCCheck = NewVar - L
        ELSE IF (NewVar.LT.0.d0) THEN
          PBCCheck = NewVar + L
        ELSE
          PBCCheck = NewVar
        ENDIF

       RETURN
      END

C-----Computes the distance for a given X and Y-------------------------
      DOUBLE PRECISION FUNCTION Module(VarX, VarY)
       IMPLICIT NONE 
       DOUBLE PRECISION VarX, VarY

       Module = DSQRT(VarX**2.d0 + VarY**2.d0)
      END

C-----Computes L -------------------------------------------------------
      DOUBLE PRECISION FUNCTION ComputesL(NPart, phi, sigma)
       IMPLICIT NONE 
       INTEGER NPart
       DOUBLE PRECISION pi, phi, sigma

       pi = ACOS(-1d0)

       ComputesL = DSQRT((pi*sigma*sigma*NPart)/(4*phi))

      END

C-----Computes MSD------------------------------------------------------
      DOUBLE PRECISION FUNCTION ComputeMSD(Xini, Yini, Xact, Yact,
     &  NPart)
       IMPLICIT NONE 
       INTEGER NPart, i
       DOUBLE PRECISION Xini(NPart), Yini(NPart), Xact(NPart),
     &  Yact(NPart), DisplacementX, DisplacementY

      ComputeMSD = 0d0

       DO i = 1, NPart
        DisplacementX = Xact(i) - Xini(i)
        DisplacementY = Yact(i) - Yini(i)
        ComputeMSD = ComputeMSD + DisplacementX*DisplacementX + 
     &    DisplacementY*DisplacementY
       ENDDO

       ComputeMSD = ComputeMSD/NPart

      END


C-----Computes Diffusivity----------------------------------------------
      DOUBLE PRECISION FUNCTION ComputeDiffusivity(MSD, t)
       IMPLICIT NONE 
       INTEGER t
       DOUBLE PRECISION MSD

      ComputeDiffusivity = MSD/(4d0*dfloat(t))

      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                            SUBROUTINES
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

C-----Computes Lx and Ly -----------------------------------------------
      SUBROUTINE ComputesLs(NPart, phi, sigma, Lx, Ly)
       IMPLICIT NONE 
       INTEGER NPart
       DOUBLE PRECISION pi, phi, sigma, Lx, Ly

       pi = ACOS(-1d0)

       Lx = DSQRT((pi*sigma*sigma*NPart)/(40*phi))
       Ly = 10*Lx
       Write(*,*) Lx, Ly

      END

C-----Initial position of NPart particles-------------------------------
      SUBROUTINE iniPositionGravity(NPart, Sigma, Lx, Ly, Xini, Yini)
      IMPLICIT NONE
      INTEGER NPart
      INTEGER i, k
      DOUBLE PRECISION Lx, Ly, Xini(NPart), Yini(NPart), 
     & DifferenceX, DifferenceY, Module, Sigma
 
      Xini(1) = Sigma/2d0 + RAND()*(Lx-Sigma)
      Yini(1) = Sigma/2d0 + RAND()*(Ly-Sigma)

      DO i = 2, NPart
20      Xini(i) = Sigma/2d0 + RAND()*(Lx-Sigma)
        Yini(i) = Sigma/2d0 + RAND()*(Ly-Sigma)

        DO k = 1, i-1

          DifferenceX = Xini(i)-Xini(k)
          DifferenceY = Yini(i)-Yini(k)

          IF (Module(DifferenceX, DifferenceY).LT.sigma) THEN
            GOTO 20
          ENDIF

        ENDDO

      ENDDO

      END

C-----Initial position of NPart particles-------------------------------
      SUBROUTINE initialPosition(NPart, Sigma, L, Xinitial, Yinitial)
      IMPLICIT NONE
      INTEGER NPart
      INTEGER i, k
      DOUBLE PRECISION L, Xinitial(NPart), Yinitial(NPart), 
     & DifferenceX, DifferenceY, Module, Sigma
 
      Xinitial(1) = RAND()*L
      Yinitial(1) = RAND()*L

      DO i = 2, NPart
21      Xinitial(i) = RAND()*L
        Yinitial(i) = RAND()*L

        DO k = 1, i-1

          DifferenceX = Xinitial(i)-Xinitial(k)
          DifferenceY = Yinitial(i)-Yinitial(k)

          IF (Module(DifferenceX, DifferenceY).LT.sigma) THEN
            GOTO 21
          ENDIF

        ENDDO

      ENDDO

      END


C-----Initial position of NPart particles-------------------------------
      SUBROUTINE initialTriangular(NPart, Sigma, L, Xini, Yini)
      IMPLICIT NONE
      INTEGER NPart
      INTEGER aux, aux2, count
      DOUBLE PRECISION L, Xini(NPart), Yini(NPart), x, y, xo, pi, Sigma

      pi = ACOS(-1d0)

      count = 1
      aux = 0

      DO WHILE (count.LE.NPart)
        y = (0.5d0 + aux*SIN(pi/3d0))*Sigma
        aux = aux + 1
        IF (MOD(aux,2).EQ.1) THEN
          xo = Sigma/2d0
        ELSE
          xo = Sigma
        ENDIF

        aux2 = 0
        x = 0d0
        DO WHILE ((x.LE.(L-Sigma*3d0/2d0)).AND.(count.LE.NPart))
          x = xo + Sigma*aux2
          Xini(count) = x
          Yini(count) = y
          count = count + 1
          aux2 = aux2 + 1
        ENDDO
      ENDDO

      END


C-----Initial position of NPart particles-------------------------------
      SUBROUTINE particleIndex(NPart, PartIndex)
      IMPLICIT NONE
      INTEGER NPart
      INTEGER i, k, AuxRand, PartIndex(NPart)
 
      DO i = 1, NPart
25      AuxRand = INT(RAND()*NPart) + 1

        DO k = 1, i-1
          IF (ANY(PartIndex == AuxRand)) THEN
            GOTO 25
          ENDIF
        ENDDO
        PartIndex(i) = AuxRand

      ENDDO

      END

C-----Fa un pas de Monte Carlo------------------------------------------
      SUBROUTINE PasMC(NPart, Sigma, L, Xini, Yini, XfinReal, YfinReal,
     &  XiniPBC, YiniPBC, XfinPBC, YfinPBC, delta)
      IMPLICIT NONE
      INTEGER i, k, NPart, PartIndex(NPart)
      DOUBLE PRECISION L, Xini(NPart), Yini(NPart), XiniPBC(NPart),
     &  Module, Sigma, delta, DiffXPBC, DiffYPBC, YiniPBC(NPart),
     &  XfinPBC(NPart), YfinPBC(NPart), XfinReal(NPart),
     &  YfinReal(NPart), IncrementX, IncrementY, PBCCheck

      XfinReal = Xini
      YfinReal = Yini

      XfinPBC = XiniPBC
      YfinPBC = YiniPBC

      CALL particleIndex(NPart, PartIndex)

      DO i = 1, NPart

        IncrementX = delta*(RAND() - 0.5d0)
        IncrementY = delta*(RAND() - 0.5d0)

        XfinReal(PartIndex(i)) = XfinReal(PartIndex(i)) + IncrementX
        YfinReal(PartIndex(i)) = YfinReal(PartIndex(i)) + IncrementY

        XfinPBC(PartIndex(i)) = PBCCheck(XfinPBC(PartIndex(i)), 
     &    IncrementX, L)
        YfinPBC(PartIndex(i)) = PBCCheck(YfinPBC(PartIndex(i)),
     &    IncrementY, L)

        DO k = 1, NPart
          IF (i.NE.k) THEN
            DiffXPBC = XfinPBC(PartIndex(i))-XfinPBC(PartIndex(k))
            DiffYPBC = YfinPBC(PartIndex(i))-YfinPBC(PartIndex(k))

            IF (Module(DiffXPBC, DiffYPBC).LT.sigma) THEN
              XfinReal(PartIndex(i)) = Xini(PartIndex(i))
              YfinReal(PartIndex(i)) = Yini(PartIndex(i))

              XfinPBC(PartIndex(i)) = XiniPBC(PartIndex(i))
              YfinPBC(PartIndex(i)) = YiniPBC(PartIndex(i))
            ENDIF
          ENDIF
        ENDDO

      ENDDO

      END


C-----Fa un pas de Monte Carlo amb walls--------------------------------
      SUBROUTINE PasMCWallsGrav(NPart, Sigma, Lx, Ly, Xini, Yini,
     &  Xfin, Yfin, delta, g, accep, total)
      IMPLICIT NONE
      INTEGER i, k, NPart, PartIndex(NPart), accep, total
      DOUBLE PRECISION Xini(NPart), Yini(NPart), Module, Sigma,
     &  delta, DiffX, DiffY, Xfin(NPart), Yfin(NPart), Lx, Ly, 
     &  IncrementX, g, IncrementY, BoltzFactor, AuxX, AuxY

      Xfin = Xini
      Yfin = Yini

      CALL particleIndex(NPart, PartIndex)

      DO i = 1, NPart

        total = total + 1
        IncrementX = delta*(RAND() - 0.5d0)
        IncrementY = delta*(RAND() - 0.5d0)

        BoltzFactor = EXP(-g*IncrementY)

        AuxX = Xfin(PartIndex(i)) + IncrementX
        AuxY = Yfin(PartIndex(i)) + IncrementY

        IF ((IncrementY.LT.0d0).OR.(BoltzFactor.GE.RAND())) THEN
        IF ((AuxX.GE.(Sigma/2d0)).AND.(AuxX.LE.(Lx-Sigma/2d0))) THEN
          IF ((AuxY.GE.(Sigma/2d0)).AND.(AuxY.LE.(Ly-Sigma/2d0))) THEN
            accep = accep + 1
            
          Xfin(PartIndex(i)) = AuxX
          Yfin(PartIndex(i)) = AuxY

          DO k = 1, NPart
          IF (i.NE.k) THEN
            DiffX = Xfin(PartIndex(i))-Xfin(PartIndex(k))
            DiffY = Yfin(PartIndex(i))-Yfin(PartIndex(k))

            IF (Module(DiffX, DiffY).LT.sigma) THEN
              Xfin(PartIndex(i)) = Xini(PartIndex(i))
              Yfin(PartIndex(i)) = Yini(PartIndex(i))
              accep = accep - 1
              EXIT
            ENDIF
          ENDIF
        ENDDO

          ENDIF
        ENDIF
        ENDIF
      ENDDO

      END

C-----Histograma--------------------------------------------------------
      SUBROUTINE Histo(NDat, XGauss, XIni, XFin, NBox, XHis, VHis, 
     &  ErrHis, IErr)

      !-----Definicio de variables input--------------------------------
      !NBOX--->Number of boxes
      !XG ---> Vector de NDAT dimensions que conté les dades per crear l'histograma
      !XIni, XFin ---> Límit inferior i límit superior respectivament

      !-----Definicio de variables output-------------------------------
      !XHIS ---> Posició central de la caixa
      !VHIS ---> Posició de la barra corresponent
      !ERRHIS ---> Estimació de l'error
      !-----------------------------------------------------------------
        IMPLICIT NONE
        INTEGER NDat,NBox,i,IBox,ICount,IErr
        DOUBLE PRECISION XGauss(NDat),XIni, XFin, BoxSize, XHis(NBox),
     &   VHis(NBox),ErrHis(NBox)
        IErr = 0
        IF (XIni.GE.XFin) THEN 
          IErr = 1
          RETURN
        ENDIF

        !Calculem el tamany que ha de tenir cada caixa
        BoxSize = (XFin - XIni)/NBox

        !Conta el nombre de punts que hi ha dins l'interval XA i XB
        ICount = 0

        !Configura totes les posicions de la caixa a 0
        DO i = 1, NBox
          VHis(i) = 0
          ErrHis(i) = 0
        ENDDO

        DO i = 1, NDat
          !Comprova si les dades estan en interval XA i XB
          IF (XGauss(i).GE.XIni.AND.XGauss(i).LE.XFin) THEN
            IBox = INT((XGauss(I) - XIni)/BoxSize) + 1
            !Posa XB a l'última caixa
            IF (IBox.EQ.NBox + 1) IBOX = NBOX
            
            VHis(IBox) = VHis(IBox) + 1
            ICount = ICount + 1
          ENDIF
        ENDDO

        IF (ICount.EQ.0) THEN 
          IErr = 2
          RETURN
        ENDIF

        IErr=0
C        PRINT*,"ACCEPTED:", ICount, " OUT OF:", NDat

        DO i = 1, NBox
          !Valor central de la caixa
          XHis(i) = XIni + BoxSize/2.D0 + (i-1)*BoxSize

          !Desviació estàndard del binomial corresponent per trobar l'estimació de l'error
          ErrHis(i) = SQRT(VHis(I)/ICount*(1.D0 - VHis(I)/ICount))/
     &    BoxSize/SQRT(DBLE(ICount))

          !Valor normalitzat de la caixa
          VHis(I) = VHis(i)/(ICount*BoxSize)
        ENDDO
       
      RETURN 
      END

C-----Exerice 1---------------------------------------------------------
      SUBROUTINE Excercise1(NPart, Sigma, L, tau)
      IMPLICIT NONE
      INTEGER NPart, i, tau, j
      DOUBLE PRECISION Xinitial(NPart), Yinitial(NPart), XfinPBC(NPart),
     & YfinPBC(NPart), XfinReal(NPart), YfinReal(NPart), MSDReal, L,
     & ComputeMSD, XiniReal(NPart), XiniPBC(NPart), YiniReal(NPart),
     & YiniPBC(NPart), Diffusivity, ComputeDiffusivity, Delta(6),
     & Sigma

      OPEN(15,FILE = "practica2_data.dat")
      OPEN(16,FILE = "practica2_data2.dat")

      Delta = (/ 0.001d0, 0.003d0, 0.01d0, 0.03d0, 0.1d0, 0.3d0 /)

      WRITE(16,*)"#   Delta                    Diffusivity"
      WRITE(15,*)"#   Time                    MDS"

      DO j = 1, SIZE(Delta)
      CALL initialPosition(NPart, Sigma, L, Xinitial, Yinitial)

      XiniReal = Xinitial
      YiniReal = Yinitial
      
      XiniPBC = Xinitial
      YiniPBC = Yinitial

      WRITE(15,*)"#  DELTA = ", Delta(j) 
      DO i = 1, tau

        CALL PasMC(NPart, Sigma, L, XiniReal, YiniReal, XfinReal, 
     &    YfinReal, XiniPBC, YiniPBC, XfinPBC, YfinPBC, Delta(j))

        MSDReal = ComputeMSD(Xinitial, Yinitial, XfinReal, YfinReal,
     & NPart)
        WRITE(15,*)i, MSDReal

        XiniReal = XfinReal
        YiniReal = YfinReal

        XiniPBC = XfinPBC
        YiniPBC = YfinPBC

      ENDDO
      WRITE(15,*)
      WRITE(15,*)

      Diffusivity = ComputeDiffusivity(MSDReal, tau)
      WRITE(16,*)Delta(j), Diffusivity

      ENDDO

      CLOSE(15)
      CLOSE(16)
      END

C-----Exerice 2---------------------------------------------------------
      SUBROUTINE Excercise2(NPart, Sigma, tau)
      IMPLICIT NONE
      INTEGER NPart, i, tau, j, k
      DOUBLE PRECISION Xinitial(NPart), Yinitial(NPart), XfinPBC(NPart),
     & YfinPBC(NPart), XfinReal(NPart), YfinReal(NPart), MSDReal, L,
     & ComputeMSD, XiniReal(NPart), XiniPBC(NPart), YiniReal(NPart),
     & YiniPBC(NPart), Delta, Sigma, Phis(3), ComputesL

      WRITE(*,*) "hola fora del bucle"
      Phis = (/ 0.05d0, 0.2d0, 0.5d0 /)
      Delta = Sigma/10d0
      OPEN(17,FILE = "practica2_E2.dat")
      OPEN(18,FILE = "practica2_E2Snapshots.dat")

      WRITE(18,*)"#           X                 Y"

      WRITE(17,*)"#      TimeStep                 MSD"
      DO k = 1, SIZE(Phis)
        WRITE(17,*)"#      Phi = ",Phis(k)

        WRITE(18,*)"#      Phi = ",Phis(k)
        L = ComputesL(NPart, Phis(k), sigma)

      CALL initialTriangular(NPart, Sigma, L, Xinitial, Yinitial)

      XiniReal = Xinitial
      YiniReal = Yinitial
      
      XiniPBC = Xinitial
      YiniPBC = Yinitial

      DO i = 0, tau

        CALL PasMC(NPart, Sigma, L, XiniReal, YiniReal, XfinReal, 
     &    YfinReal, XiniPBC, YiniPBC, XfinPBC, YfinPBC, Delta)

        MSDReal = ComputeMSD(Xinitial, Yinitial, XfinReal, YfinReal,
     & NPart)
        WRITE(17,*)i, MSDReal

        XiniReal = XfinReal
        YiniReal = YfinReal

        XiniPBC = XfinPBC
        YiniPBC = YfinPBC

        IF (MOD(i, 20000).EQ.0) THEN
          WRITE(18,*)"#temps = ", i
          WRITE(*,*) "hola"
          DO j = 1, NPart
            WRITE(18,*)XfinPBC(j), YfinPBC(j)
          ENDDO
          WRITE(18,*)
          WRITE(18,*)
        ENDIF

      ENDDO
      WRITE(*,*) "hola"

      WRITE(17,*)
      WRITE(17,*)
        ENDDO
      CLOSE(17)
      CLOSE(18)

      END

C-----Exerice 3, 4 and 5------------------------------------------------
      SUBROUTINE Excercise5(NPart, Sigma, phi, tau)
      IMPLICIT NONE
      INTEGER NPart, i, tau, j, k, accep, total, box, m, ttau
      PARAMETER(box = 20)
      DOUBLE PRECISION Xinitial(NPart), Yinitial(NPart), 
     & Xfin(NPart), Yfin(NPart), MSD, Lx, Ly, ComputeMSD, Xini(NPart),
     & Yini(NPart), Delta, Sigma, Gravity(5), phi, h, rho(box)

      OPEN(35,FILE = "practica2_E5.dat")
      OPEN(36,FILE = "practica2_E5Snap.dat")
      OPEN(14, FILE = "rhos.dat")

      Delta = 0.5d0
      Gravity = (/ 0d0, 0.01d0, 0.1d0, 1d0, 10d0 /)

      WRITE(36,*)"#      X                     Y"
      WRITE(35,*)"#    Time                    MDS"

      CALL ComputesLs(NPart, phi, Sigma, Lx, Ly)

      h = Ly/box

      DO k = 1, SIZE(Gravity)
        accep = 0
        total = 0

        WRITE(36,*)"#      Gravity = ", Gravity(k)

        WRITE(35,*)"#      Gravity = ", Gravity(k)
      CALL iniPositionGravity(NPart, Sigma, Lx, Ly, Xinitial, Yinitial)

      Xini = Xinitial
      Yini = Yinitial
      WRITE(*,*)"hola"

      IF (DABS(Gravity(k)-10d0).LE.10.**(-5.)) THEN
        ttau = 1000
      ELSE
        ttau = 5000
      ENDIF
      WRITE(*,*)ttau

      DO i = 1, ttau

        CALL PasMCWallsGrav(NPart, Sigma, Lx, Ly, Xini, Yini, 
     &    Xfin, Yfin, Delta, Gravity(k), accep, total)

        MSD = ComputeMSD(Xinitial, Yinitial, Xfin, Yfin,
     &    NPart)
        WRITE(35,*)i, MSD

        Xini = Xfin
        Yini = Yfin

      ENDDO

      WRITE(35,*)
      WRITE(35,*)

      DO m = 1, box
        rho(m) = 0.d0
      ENDDO

      DO j = 1, NPart

        DO m = 1, box
          IF (Yfin(j).LT.(m*h)) THEN
            rho(m) = rho(m) + 1.d0/DBLE(NPart)
            EXIT
          ENDIF
      ENDDO

        WRITE(36,*)Xfin(j), Yfin(j)
      ENDDO

      WRITE(36,*)
      WRITE(36,*)
      WRITE(*,*)"# Gavity = ", Gravity(k), "Trans prob = ", 
     &dble(accep)/dble(total), accep

      DO m = 1, box
        WRITE(14,*) m*h, rho(m)
      ENDDO
      WRITE(14,*)
      WRITE(14,*)
      
      ENDDO
       
      CLOSE(35)
      CLOSE(36)
      CLOSE(14)
      END
