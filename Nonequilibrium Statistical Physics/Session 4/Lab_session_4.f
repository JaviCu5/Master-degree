C23456789012345678901234567890123456789012345678901234567890123456789012
c-----PROJECT 3: Continuous-Time Monte Carlo or Gillespie Algorithm-----
c-----Javier Castillo Uviña
      PROGRAM Practica3
      IMPLICIT NONE
      INTEGER seed, NPart, TotalStep, vecTau(3), Ndat, i, q, k, h,
     &  realizations, taus(5)
      DOUBLE PRECISION DeltaT, L
      PARAMETER (NPart = 200, TotalStep = 1000, DeltaT = 0.1,
     &  NDat = 5*10**6)
      DOUBLE PRECISION  XGauss(NDat), Diffus, vecMSD(TotalStep),
     &  vecDisp(TotalStep), forces(5), D, mu, force

      seed = 1
      Call srand(seed)
      vecTau = (/1, 10, 100/)

      L = 20d0
      Diffus = 3.d0
c-----------------------------------------------------------------------
      !SECTION 1
      !Random walk with memory
c-----------------------------------------------------------------------        

      OPEN(14, FILE = "MSD.dat")
      OPEN(15, FILE = "Autocorr.dat")
      OPEN(19, FILE = "Trajectories.dat")

      CALL RandGauss(NDat, 0d0, 1d0, XGauss)
      DO i = 1, SIZE(vecTau)
        CALL RandWalk(Npart, XGauss, NDat, L, TotalStep, DeltaT,
     &    vecTau(i), Diffus, 14, 15)
        WRITE(14,*)
        WRITE(14,*)
        WRITE(15,*)
        WRITE(15,*)
        WRITE(19,*)
        WRITE(19,*)
      ENDDO

      CLOSE(14)
      CLOSE(15)
      CLOSE(19)

c-----------------------------------------------------------------------
      !SECTION 2
      !Correlation function
c----------------------------------------------------------------------- 

      forces = (/ 0., 0.001, 0.01, 0.1, 1. /)
      realizations = 1

      OPEN(16, FILE = "XDisp.dat")

      DO k = 1, SIZE(forces)
        DO q = 1, TotalStep
          vecMSD(q) = 0d0
          vecDisp(q) = 0d0
        ENDDO

        DO h = 1, realizations
          seed = 4 + seed
          CALL srand(seed)
          CALL RandGauss(NDat, 0d0, 1d0, XGauss)

          CALL ForceRandWalk(Npart, XGauss, NDat, L, TotalStep, DeltaT, 
     &      vecTau(2), Diffus, forces(k), vecMSD, vecDisp, realizations)
          vecMSD = vecMSD
          vecDisp = vecDisp
        ENDDO

        DO q = 1, TotalStep
          WRITE(16,*)q*DeltaT, vecMSD(q), vecDisp(q)
        ENDDO

        WRITE(16,*)
        WRITE(16,*)

      ENDDO

      CLOSE(16)

c-----------------------------------------------------------------------
      !SECTION 3
c----------------------------------------------------------------------- 

      OPEN(17, FILE = "muDiff.dat")

      force = forces(5)
      taus = (/1, 5, 10, 50, 100/)

      DO k = 1, SIZE(taus)
        DO q = 1, TotalStep
          vecMSD(q) = 0d0
          vecDisp(q) = 0d0
        ENDDO

        DO h = 1, realizations
          seed = 4 + seed
          CALL srand(seed)
          CALL RandGauss(NDat, 0d0, 1d0, XGauss)

          
          CALL ForceRandWalk(Npart, XGauss, NDat, L, TotalStep, DeltaT, 
     &      taus(k), Diffus, force, vecMSD, vecDisp, realizations)
          vecMSD = vecMSD
          vecDisp = vecDisp
        ENDDO

        mu = vecDisp(TotalStep-1)/(TotalStep*DeltaT*force)
        D = vecMSD(TotalStep-1)/(TotalStep*DeltaT*2.d0)
        WRITE(17,*)taus(k), mu, D
      ENDDO

      CLOSE(17)

      Call system ("gnuplot -p script.gnu")
      END PROGRAM

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                             FUNCTIONS
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

C-----RandWalk step-----------------------------------------------------
      DOUBLE PRECISION FUNCTION VarRand(Var, Incr)
      IMPLICIT NONE
      DOUBLE PRECISION Var, Incr

      VarRand = Var + Incr
      END

C-----RandWalk step with PBC--------------------------------------------
      DOUBLE PRECISION FUNCTION VarPBC(L, Var, Incr)
      IMPLICIT NONE
      DOUBLE PRECISION Var, Incr, L, NewVar

      NewVar = Var + Incr
      IF (NewVar.GT.L) THEN
        VarPBC = NewVar - L
      ELSE IF (NewVar.LT.0.d0) THEN
        VarPBC = NewVar + L
      ELSE
        VarPBC = NewVar
      ENDIF

      RETURN
      END

C-----Computes the distance for a given X and Y-------------------------
      DOUBLE PRECISION FUNCTION Modul(VarX, VarY)
      IMPLICIT NONE 
      DOUBLE PRECISION VarX, VarY

      Modul = DSQRT(VarX**2.d0 + VarY**2.d0)
      END

C-----Computes the distance for a given X and Y-------------------------
      DOUBLE PRECISION FUNCTION MSDAnalitic(D, t, tau)
      IMPLICIT NONE 
      DOUBLE PRECISION D, t, tau

      MSDAnalitic = 4*D*(t + tau*(DEXP(-t/tau) - 1))
      END

C-----Computes the distance for a given X and Y-------------------------
      DOUBLE PRECISION FUNCTION AutocorrAnalitic(D, t, tau)
      IMPLICIT NONE 
      DOUBLE PRECISION D, t, tau

      AutocorrAnalitic = (D/tau) * DEXP(-t/tau)
      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                            SUBROUTINES
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


C-----Random gaussian N number generator with Mu and Sigma parameters---
      SUBROUTINE RandGauss (NDat, Mu, Sigma, XGauss)
      IMPLICIT NONE
      DOUBLE PRECISION Rand1, Rand2, R, Pi, Phi, XGauss(NDat), X1, X2,
     &  Sigma, Mu
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

C-----MSD for N particles in a time T=DeltaT*TotalSt for a specific tau-
      SUBROUTINE RandWalk(Npart, XGauss, NDat, L, TotalStep, DeltaT,
     &  VTau, D, iFile, iFIle2)
      IMPLICIT NONE
      INTEGER Ndat, i, j, NPart, TotalStep, iFile, VTau, iFIle2
      DOUBLE PRECISION XIo(NPart), YIo(NPart), XI(NPart), YI(NPart),
     &  XGauss(NDat), IncrX, IncrY, MSD, VarRand, VarPBC, L,
     &  XIPBC(NPart), YIPBC(NPart), D, Modul, tau, DeltaT, MSDAna,
     &  MSDAnalitic, XiX(NPart), XiY(NPart), XiXo(NPart), XiYo(NPart),
     &  AutocorrAnalitic, Autocorr, AutocorrAn

      tau = DBLE(VTau)

      DO j = 1, Npart
        XIo(j) = RAND()*DBLE(L)
        YIo(j) = RAND()*DBLE(L)

        XI(j) = XIo(j)
        YI(j) = YIo(j)

        XiX(j) = DSQRT(D/tau)
        XiY(j) = DSQRT(D/tau)

        XiXo(j) = XiX(j)
        XiYo(j) = XiY(j)
      ENDDO

      DO i = 1, TotalStep

        MSD = 0d0
        MSDAna = MSDAnalitic(D, i*DeltaT, tau)
        Autocorr = 0d0
        AutocorrAn = AutocorrAnalitic(D, i*DeltaT, tau)
      
        DO j = 1, NPart
          XiX(j) = XiX(j)*(1-DeltaT/tau)+DSQRT((2*D*DeltaT)/tau**2)*
     &      XGauss(2*j + i*TotalStep)
          IncrX = XiX(j)*DeltaT

          XiY(j) = XiY(j)*(1-DeltaT/tau)+DSQRT((2*D*DeltaT)/tau**2)*
     &      XGauss(2*j + 1 + i*TotalStep)
          IncrY = XiY(j)*DeltaT

          XI(j) = VarRand(XI(j), IncrX)
          YI(j) = VarRand(YI(j), IncrY)

          XIPBC(j) = VarPBC(L, XIPBC(j), IncrX)
          YIPBC(j) = VarPBC(L, YIPBC(j), IncrY)

          MSD = MSD + Modul(XI(j)-XIo(j), YI(j)-YIo(j))**2d0

          Autocorr = Autocorr + XiXo(j)*XiX(j) + XiYo(j)*XiY(j)
          IF (j.EQ.1) THEN
            WRITE(19,*)XI(j), YI(j)
          ENDIF
        ENDDO

        MSD = MSD/DBLE(NPart)
        Autocorr = Autocorr/DBLE(NPart*2)

        WRITE(iFile,*) i*DeltaT, MSD, MSDAna
        WRITE(iFIle2,*) i*DeltaT, Autocorr, AutocorrAn

      ENDDO

      RETURN
      END


C-----MSD for N particles in a time T=DeltaT*TotalSt for a specific tau-
C-----in presence of a force--------------------------------------------
      SUBROUTINE ForceRandWalk(Npart, XGauss, NDat, L, TotalStep, 
     &  DeltaT, VTau, D, f, vecMSD, vecDisp, realisations)
      IMPLICIT NONE
      INTEGER Ndat, i, j, NPart, TotalStep, VTau, realisations
      DOUBLE PRECISION XIo(NPart), XI(NPart), XiX(NPart), f,
     &  XGauss(NDat), IncrX, MSD, VarRand, L, coefF(NPart), disp, 
     &  D, tau, DeltaT, vecMSD(TotalStep), vecDisp(TotalStep)

      tau = DBLE(VTau)

      DO j = 1, Npart
        XIo(j) = RAND()*DBLE(L)

        XI(j) = XIo(j)

        XiX(j) = DSQRT(D/tau)

        IF (RAND().LT.0.5d0) THEN
            coefF(j) = -1d0
        ELSE
            coefF(j) = 1d0 
        ENDIF
      ENDDO

      DO i = 1, TotalStep
        MSD = 0d0
        disp = 0d0
      
        DO j = 1, NPart
          XiX(j) = XiX(j)*(1-DeltaT/tau)+DSQRT((2*D*DeltaT)/(tau**2.))*
     &      XGauss(j + i*TotalStep)
          IncrX = (coefF(j)*f + XiX(j))*DeltaT

          XI(j) = VarRand(XI(j), IncrX)

          MSD = MSD + (XI(j)-XIo(j))**2d0
          disp = disp + coefF(j)*(XI(j)-XIo(j))

        ENDDO

        vecMSD(i) = vecMSD(i) + MSD/DBLE(NPart*realisations)
        vecDisp(i) = vecDisp(i) + disp/DBLE(NPart*realisations)

      ENDDO
      WRITE(*,*)vecDisp(TotalStep-1)

      RETURN
      END