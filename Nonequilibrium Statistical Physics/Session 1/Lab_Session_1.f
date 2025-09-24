c-----PRÀCTICA 1: Brownian Dynamics, Fluctuations and Response----------

c-----Javier Castillo Uviña


c-----MAIN PROGRAM------------------------------------------------------
      PROGRAM project1


c-----Variables---------------------------------------------------------
      IMPLICIT NONE
      !SECTION 2
      INTEGER ISeed, NDat, i, NPart, TimeSteps, VecT(5), IErr
      DOUBLE PRECISION Sigma, Mu
      PARAMETER (NDat = 10**8)
      PARAMETER (NPart = 1000)
      PARAMETER (TimeSteps = 100)
      PARAMETER (ISeed = 16866231)
      DOUBLE PRECISION XGauss(NDat)
      INTEGER NBox
      PARAMETER (NBox = 100)
      DOUBLE PRECISION A, B, XHis(nbox), VHis(nbox), ErrHis(nbox),
     &  DeltaT, Sum1, Sum2, Desv, X, Y, XR , YR, L, tau,
     &  Diffusivity, Gammas(5)
      PARAMETER (DeltaT = 0.1d0)
      PARAMETER (L = 100.d0)
      PARAMETER (tau = TimeSteps*DeltaT)
      Gammas = (/ 0.1d0, 1d0, 3d0, 10d0, 30d0 /)

      Mu = 0.d0
      Sigma = 1.d0

c-----------------------------------------------------------------------
      !SECTION 2
c-----------------------------------------------------------------------
      !Generate ranDOm numbers n following a normal distribution
      CALL SRand(ISeed)
      CALL RandGauss (NDat, Mu, Sigma, XGauss)

      !Distribution
      A = - 4 * Sigma
      B = 4 * Sigma
      CALL Histo(NDat, XGauss, A, B, NBox, XHis, VHis, ErrHis, IErr)


      OPEN(15,FILE = "project1-res.dat")
      WRITE(15,*) "# Normal distribution"
      WRITE(15,'(a20,2x,a20,4x,a20)') "# XHIS(I)", "VHIS(I)", "Error"
      DO i = 1, NBox
        WRITE(15,*) XHis(i), VHis(i), ErrHis(i)
      ENDDO
      
      WRITE(15,*)
      WRITE(15,*)


      Sum1=0.d0
      Sum2=0.d0

      DO i = 1, Ndat
        X = XGauss(i)
        Sum1 = Sum1 + x
        Sum2 = Sum2 + x*x
      ENDDO

      !Mitjana
      Sum1 = Sum1/REAL(NDat)
      Sum2 = Sum2/REAL(NDat)

      !Desviació
      Desv = DSQRT(Sum2 - Sum1*Sum1)

      WRITE(*,*) "Mitjana", Sum1, "Desviació", Desv
      !Correct results mitjana = 0 AND sigma = 1


c-----------------------------------------------------------------------
      !SECTION 3
c-----------------------------------------------------------------------
      WRITE(15,*) "# X Coordinate","  Y Coordinate"
      XR = 0.d0
      YR = 0.d0
      WRITE(15,*) XR, YR
      
      DO i = 1, TimeSteps
        X = XGauss(2*i-1)
        Y = XGauss(2*i)
        XR = XR + X*DSQRT(2*Gammas(3)*DeltaT)
        YR = YR + Y*DSQRT(2*Gammas(3)*DeltaT)
        WRITE(15,*) XR, YR
      ENDDO

      WRITE(15,*)
      WRITE(15,*)
      

c-----------------------------------------------------------------------
      !SECTION 4 and 5
c-----------------------------------------------------------------------


      DO i = 1, SIZE(Gammas)

      CALL RandWalkPart(Ndat, Npart, XGauss, L, TimeSteps, DeltaT,
     &  Gammas(i), tau, Diffusivity)

        IF (i.EQ.1) THEN
          WRITE(15,*) "#Diffusivity         Gamma"
        ENDIF
        WRITE(15,*) Diffusivity, Gammas (i)

      ENDDO
      
      CLOSE(15)

c-----------------------------------------------------------------------
      !SECTION 6
c-----------------------------------------------------------------------

      VecT = (/ 1, 10, 30, 100, 1000 /)

      CALL DiffEquation(Ndat, NPart, XGauss, L, DeltaT,
     &  Gammas(3), VecT)


      CALL system ("gnuplot -p Graphics.gnu")

      END PROGRAM project1
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
      DOUBLE PRECISION FUNCTION CalculateR(VarX, VarY)
       IMPLICIT NONE 
       DOUBLE PRECISION VarX, VarY

       CalculateR = DSQRT(VarX**2.d0 + VarY**2.d0)
      END


C-----Random gaussian number generator with Mu and Sigma parameters
      SUBROUTINE RandGauss (NDat, Mu, Sigma, XGauss)
        IMPLICIT NONE
        DOUBLE PRECISION Rand1, Rand2, R, Pi, Phi, XGauss(NDat), X1, X2,
     &   Sigma, Mu
        INTEGER NDat,i

        Pi = ACOS(-1.d0)
        
        !Box-Muller algorithm
        DO i = 1, NDat - 1, 2
        !With rand() function we generate 2 ranDOm nombers between 0 and 1
              Rand1 = RAND()
              Rand2 = RAND()

            !Variables changes to UnIForm-Gaussian distribution
              R = DSQRT(-2.d0*log(Rand1))
              Phi = 2*Pi*Rand2

            !Obtain two values that follow a general gaussian
              X1 = R*COS(Phi)
              X2 = R*SIN(Phi)

            !Obtain xmuller distributed ith the gaussian we want
            !If mu =! 0 and sigma!=1 we have to DO a variable change 
              XGauss(i) = Mu + X1*Sigma
              XGauss(i+1) = Mu + X2*Sigma
        ENDDO
      RETURN
      END


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                            SUBROUTINES
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


C-----We will obtain the particles distribution-------------------------
      SUBROUTINE DiffEquation(Ndat, NPart, XGauss, L, DeltaT,
     &  Gamma, VecT)
      IMPLICIT NONE
      INTEGER Ndat, TimeSteps, i, j, NPart, k, NBox, VecT(5), IErr,
     & Count, q
      PARAMETER(NBox = 50)
      DOUBLE PRECISION XFPo(NPart), L, XGauss(NDat), XFPI(NPart),
     &   X, PBCCheck, DeltaT, Gamma, DeltaX(NPart),XHis(NBox),
     &   VHis(NBox), ErrHis(NBox)

      q = 1
      TimeSteps = INT(VecT(5)/DeltaT)
      WRITE(*,*) "TimeSteps ", TimeSteps

      OPEN(16,FILE = "FPHistograms.dat")

      DO i = 1, NPart
        XFPo(i) = RAND()*L
      ENDDO

      DO i = 1, NPart
          X = XGauss(i)
          XFPI(i) = PBCCheck(XFPo(i), X*DSQRT(2*Gamma*DeltaT), L)
          DeltaX(i) = XFPI(i) - XFPo(i)
        ENDDO

      DO j = 2, TimeSteps

        Count = INT(j*DeltaT)
        
        IF (VecT(q).EQ.count) THEN
          q = q + 1
          CALL Histo(NPart, DeltaX, -L, L, NBox, XHis, VHis, 
     &    ErrHis, IErr)
          WRITE(16,*)"# Ditribution for t = ", Count
          WRITE(16,'(a20,2x,a20,4x,a20)') "# XHIS(I)", "VHIS(I)",
     &    "Error"
          DO k = 1, NBox
            WRITE(16,*) XHis(k), VHis(k), ErrHis(k)
          ENDDO
          WRITE(16,*)
          WRITE(16,*)
        ENDIF

        DO i = 1, NPart
          X = XGauss(i + (j-1)*TimeSteps)
          XFPI(i) = PBCCheck(XFPI(i), X*DSQRT(2*Gamma*DeltaT), L)
          DeltaX(i) = XFPI(i) - XFPo(i)
        ENDDO

      ENDDO

      CLOSE(16)
      RETURN
      END


      SUBROUTINE RandWalkPart(Ndat, Npart, XGauss, L, TimeSteps, DeltaT,
     &  Gamma, tau, Diffusivity)
      IMPLICIT NONE
      INTEGER Ndat, TimeSteps, i, j, NPart
      DOUBLE PRECISION XIo, YIo, L, RIo, RI, XGauss(NDat), XI, YI,
     &   NetDisplacement, X, Y, PBCCheck, DeltaT, CalculateR, Gamma,
     &   tau, NetDisplacementAverage, Diffusivity, t

      NetDisplacementAverage = 0.d0

      DO i = 1, NPart
        XIo = RAND()*L
        YIo = RAND()*L

        RIo = CalculateR(XIo, YIo)
        NetDisplacement = 0.d0

        X = XGauss(2-1 + (i-1)*2*TimeSteps)
        Y = XGauss(2 + (i-1)*2*TimeSteps)
        XI = PBCCheck(XIo, X*DSQRT(2*Gamma*DeltaT), L)
        YI = PBCCheck(YIo, Y*DSQRT(2*Gamma*DeltaT), L)

        RI = CalculateR(XI, YI)
        NetDisplacement = NetDisplacement + (RI-RIo)**2

        DO j = 2, TimeSteps
          X = XGauss(2*j-1 + (i-1)*2*TimeSteps)
          Y = XGauss(2*j + (i-1)*2*TimeSteps)
          XI = PBCCheck(XI, X*DSQRT(2*Gamma*DeltaT), L)
          YI = PBCCheck(YI, Y*DSQRT(2*Gamma*DeltaT), L)

          RI = CalculateR(XI, YI)
          NetDisplacement = NetDisplacement + (RI-RIo)**2
        ENDDO
        NetDisplacement = NetDisplacement/TimeSteps

        IF (Gamma.EQ.(0.1d0)) THEN
          IF (i.EQ.1) THEN
          WRITE(15,*) "# This data will be the first and the last",
     &    " coordinates for each particle, separated for a blank line"
          WRITE(15,*) "# X0 Coordinate","     Y0 Coordinate",
     &    "       Xf Coordinate            Yf Coordinate"
          ENDIF

          WRITE(15,*) XIo, YIo, XI, YI
          WRITE(15,*)
          IF (i.EQ.NPart) THEN
          WRITE(15,*)
          ENDIF

        ENDIF
      

      NetDisplacementAverage = NetDisplacementAverage + 
     &   NetDisplacement

      ENDDO

      NetDisplacementAverage = NetDisplacementAverage/NPart

      t = tau
      Diffusivity = NetDisplacementAverage/(4*t)
      RETURN
      END



      SUBROUTINE Histo(NDat, XGauss, XIni, XFin, NBox, XHis, VHis, 
     &  ErrHis, IErr)

      !-----Definicio de variables input-----------------------------------------------------------------------
      !NBOX--->Number of boxes
      !XG ---> Vector de NDAT dimensions que conté les dades per crear l'histograma
      !XIni, XFin ---> Límit inferior i límit superior respectivament

      !-----Definicio de variables output-----------------------------------------------------------------------
      !XHIS ---> Posició central de la caixa
      !VHIS ---> Posició de la barra corresponent
      !ERRHIS ---> Estimació de l'error
      !--------------------------------------------------------------------------------------------------------
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
