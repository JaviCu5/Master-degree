C23456789012345678901234567890123456789012345678901234567890123456789012
c-----PROJECT 3: Continuous-Time Monte Carlo or Gillespie Algorithm-----
c-----Javier Castillo Uviña
      PROGRAM Practica3
      IMPLICIT NONE
      INTEGER seed, L1, i, j, L, N, MCS, varJ, iSeed, u
      PARAMETER (L1 = 50, L = 10, N =L*L, MCS = 1000)
      INTEGER S1(1:L1, 1:L1), S(1:L, 1:L), PBC(0:L+1), vecEps(5),
     &  C(1:L,1:L), degen(5), AuxS(1:L, 1:L)
      DOUBLE PRECISION vecMagFMC(MCS), vecTransProb(5), vecMagCTMC(MCS),
     &  meanMagFMC(MCS), meanMagCTMC(MCS)

      seed = 4
      Call srand(seed)
      varJ = 1
      vecEps = (/-8*varJ, 8*varJ, -4*varJ, 4*varJ, 0/)

      !PBC connects the first and the last spin of each row and column
      PBC(0) = L
      PBC(L + 1) = 1
      DO i = 1, L
        PBC(i) = i
      ENDDO
c-----------------------------------------------------------------------
      !SECTION 1
      !Construct a disordered configuration of the Ising model
c-----------------------------------------------------------------------        

      !Create the matrix
      OPEN(14, FILE="IniMatrix.dat")
      CALL InitialConfiguration(L1, S1)

      DO i = 1, L1
        DO j = 1, L1
          WRITE(14,*)i, j, S1(i, j)
        ENDDO
        WRITE(14,*)
      ENDDO
      CLOSE(14)

      !Calculate the probabilities
      Call TransitionProb(vecEps, L, vecTransProb)

      OPEN(15, FILE="FSMCMagne.dat")
      OPEN(16, FILE="CTMCMagne.dat")

c-----------------------------------------------------------------------
      !SECTION 4
      !Do it for 10 seeds
c-----------------------------------------------------------------------    
      
      DO u = 1, MCS
        meanMagFMC(u) = 0.d0
        meanMagCTMC(u) = 0.d0
      ENDDO

      DO iSeed = 1, 10
      seed = seed + iSeed
      Call srand(seed)

c-----------------------------------------------------------------------
      !SECTION 2
      !Implement the standard Metropolis Monte Carlo scheme
c-----------------------------------------------------------------------      
      Call InitialConfiguration(L, S)
      AuxS = S

      WRITE(15,*)'# Metropolis FSMC Seed =', seed
      CALL MetropolisFSMC(L, N, MCS, S, PBC, 15, vecMagFMC)
      WRITE(15,*)
      WRITE(15,*)

c-----------------------------------------------------------------------
      !SECTION 3
      !Implement the CTMC for the same FM Ising model with PBC 
c-----------------------------------------------------------------------  
      !Classification and heatMap 
      Call ClassificationMatrix(L, S, C, PBC, degen, vecEps)

      WRITE(16,*)'# MetropolisCTMC Seed =', seed
      Call MetropolisCTMC(L, N, MCS, PBC, vecEps, AuxS, vecTransProb,
     & vecMagCTMC, 16)

      WRITE(16,*)
      WRITE(16,*)

      DO u = 1, MCS
        meanMagFMC(u) = meanMagFMC(u) + vecMagFMC(u)
        meanMagCTMC(u) = meanMagCTMC(u) + vecMagCTMC(u)
      ENDDO

      ENDDO

      WRITE(15,*)'# Mean FSMC algorithm'
      WRITE(16,*)'# Mean CTMC algorithm'

      DO u = 1, MCS
        WRITE(15,*)u, meanMagFMC(u)/DBLE(10)
        WRITE(16,*)u, meanMagCTMC(u)/DBLE(10)
      ENDDO

      CLOSE(15)
      CLOSE(16)

      Call system ("gnuplot -p s.gnu")
      END PROGRAM
      
      


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                             FUNCTIONS
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

C-----Gets the energy for a specific location of a spin matrix----------
      INTEGER FUNCTION LocEnergy(S, L, i, j, PBC)
      !-----Variables---------------------------------------------------
            !S -> Spin matrix
            !L -> Matrix lateral length
            !i -> Vertical index
            !j -> Horizontal index
            !PBC -> Periodic boundary condition vector
      !-----------------------------------------------------------------
      IMPLICIT NONE    
      INTEGER S(1:L, 1:L), i, j, L, PBC(0:L+1)

      LocEnergy = 2*S(i,j)*(S(i,PBC(j-1)) + S(i,PBC(j+1)) +
     &  S(PBC(i-1),j) + S(PBC(i+1),j))
      END

C-----Gets the spin matrix magnetization--------------------------------
      DOUBLE PRECISION FUNCTION magne(S, L, N)
      !-----Variables---------------------------------------------------
            !S -> Spin matrix
            !L -> Matrix lateral length
            !N -> LxL value
      !-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER S(1:L, 1:L), i, j, L, N
      magne = 0.D0
      DO i = 1, L
        DO j = 1, L
            magne = magne + S(i,j)
        ENDDO
      ENDDO
      magne = magne/dble(N)
      END

C-----Gets the index of the N disks randomly----------------------------
      INTEGER FUNCTION spinClass(i, j, L, S, vecEps, PBC)
      !-----Variables---------------------------------------------------
            !i -> Vertical index
            !j -> Horizontal index
            !L -> Matrix lateral length      
            !S -> Spin matrix
            !vecEps -> Are the diff spin configurations energy contribution
            !PBC -> Periodic boundary condition vector
      !-----------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER i, j, L, S(1:L, 1:L), PBC(0:L+1), LocEnergy, energy, 
     &  vecEps(5)
      
      energy = LocEnergy(S, L, i, j, PBC)
      IF (energy.EQ.(vecEps(1))) THEN
            spinClass = 1
      ELSE IF (energy.EQ.(vecEps(2))) THEN
            spinClass = 2
      ELSE IF (energy.EQ.(vecEps(3))) THEN
            spinClass = 3
      ELSE IF (energy.EQ.(vecEps(4))) THEN
            spinClass = 4
      ELSE
            spinClass = 5
      ENDIF
      END

C-----Chose a l spin-type to flip---------------------------------------
      INTEGER FUNCTION ClassChoice(vecTransProb, degen)
      !-----Variables---------------------------------------------------
      !vecTransProb -> Transition probabilities for each type
      !degen -> Vector of degeneration for each type of transition
      !-----------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER i, degen(5)
      DOUBLE PRECISION vecTransProb(5), vecProb(5), suma, random
      suma = 0.d0
      DO i = 1, 5
        suma = suma + degen(i)*vecTransProb(i)
      ENDDO

      DO i = 1, 5
        vecProb(i) = (degen(i)*vecTransProb(i))/suma
      ENDDO

      random = RAND()
      IF (random.LE.vecProb(1)) THEN
        Classchoice = 1
      ELSEIF (random.LE.vecProb(3)) THEN
        Classchoice = 3
      ELSE
        ClassChoice = 5
      ENDIF
      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c                            SUBROUTINES
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

C-----Get de initial Spin values for i and j position-------------------
      SUBROUTINE InitialConfiguration(L, S)
      !-----Variables---------------------------------------------------
            ! INPUTS
            !L -> S lenght
            ! OUTPUTS
            !S -> Spin Matrix
      !-----------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER i, j, L, S(1:L, 1:L)
      DO i = 1, L
        DO j = 1, L
          IF (RAND().lt.0.5D0) THEN
            S(i,j) = 1
          ELSE
            S(i,j) = -1
          ENDIF
        ENDDO
      ENDDO
      END

C-----Computes the magnetization evolution for FSMS---------------------
      SUBROUTINE MetropolisFSMC(L, N, MCS, S, PBC, iFile, vecMag)
      !-----Variables---------------------------------------------------
            ! INPUTS
            !L -> S lateral lenght
            !N -> LxL Spins
            !MCS -> Total MonteCarlo steps
            !S -> Spin Matrix
            !PBC -> Periodic Boundary Conditions vector
            !iFile -> The file number where it writes each magnetization
            ! OUTPUTS
            !vecMag -> Saves S magnetization for each MC step
      !-----------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER i, j, iPas, iMC, L, N, MCS, LocEnergy, IncEnergy, 
     &  S(1:L,1:L), PBC(0:L+1), iFile
      DOUBLE PRECISION mag, magne, vecMag(MCS)

      iMC = 1
      mag = DABS(magne(S, L, N))
      WRITE(iFile,*)iMC, mag
      vecMag(iMC) = mag

      DO iMC = iMC + 1, MCS

        DO iPas = 1, N 
          i = INT(L*RAND()) + 1
          j = INT(L*RAND()) + 1
          !Computes the energy for a spin change
          IncEnergy = LocEnergy(S, L, i, j, PBC)
          !Energy goes down so we accept the change
          IF (IncEnergy.le.0) THEN
                S(i,j) = -S(i,j)
          ENDIF
        ENDDO

        mag = DABS(magne(S, L, N))
        WRITE(iFile,*)iMC, mag
        vecMag(iMC) = mag

      ENDDO

      RETURN
      END

C-----Gets the transition probabilities-----------------------------
      SUBROUTINE TransitionProb(vecEps, L, vecTransProb)
      !-----Variables---------------------------------------------------
            ! INPUTS
            !vecEps -> Are the diff spin configurations energy contribution
            !L -> Matrix lateral lenght
            ! OUTPUTS
            !vecTransProb -> Transition probabilities for each type
      !-----------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER i, L, vecEps(5)
      DOUBLE PRECISION vecTransProb(5)
      DO i = 1, SIZE(vecEps)
        IF (vecEps(i).LE.0)THEN
          vecTransProb(i) = 1/DBLE(L)**2.d0
        ELSE
          vecTransProb(i) = 0.D0
        ENDIF
      ENDDO

      RETURN
      END

C-----Computes a type-spin matrix and a degeneration transition-type vector
      SUBROUTINE ClassificationMatrix(L, S, C, PBC, degen, vecEps)
      !-----Variables---------------------------------------------------
            ! INPUTS
            !L -> Matrix lateral lenght
            !S -> Spin Matrix
            !PBC -> Periodic Boundary Conditions vector
            !vecEps -> Are the diff spin configurations energy contribution
            ! OUTPUTS
            !degen -> vector of degeneration for each type of transition
            !C -> Classification matrix
      !-----------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER i, j, L, S(1:L, 1:L), C(1:L, 1:L), PBC(0:L+1), spinType,
     &  spinClass, degen(5), vecEps(5)

      degen = (/0, 0, 0, 0, 0/)

      DO i = 1, L
        DO j = 1, L
          spinType = spinClass(i, j, L, S, vecEps, PBC)
          C(i, j) = vecEps(spinType)
          degen(spinType) = degen(spinType) + 1
        ENDDO
      ENDDO

      RETURN
      END

C-----Computes the magnetization evolution for CTMC---------------------
      SUBROUTINE MetropolisCTMC(L, N, MCS, PBC, vecEps, S, vecTransProb,
     &  vecMag, iFile)
      !-----Variables---------------------------------------------------
            ! INPUTS
            !L -> Matrix lateral lenght
            !N -> LxL Spins
            !MCS -> Total MonteCarlo steps
            !PBC -> Periodic Boundary Conditions vector
            !vecEps -> Are the diff spin configurations energy contribution
            !S -> Spin Matrix
            !vecTransProb -> Transition probabilities for each type
            !iFile -> The file number where it writes
            ! OUTPUTS
            !vecMag -> Saves S magnetization for each MC step
      !-----------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER L, N, ClassChoice, spinType, degen(5), g_l, MCS, spinK,
     &  MC, S(1:L, 1:L), C(1:L, 1:L), PBC(0:L+1), vecEps(5), Delta_t,
     &  iMax, i, iFile
      DOUBLE PRECISION magne, mag, vecTransProb(5), Q, vecMag(MCS)
 
      MC = 1
      mag = DABS(magne(S, L, N))
      vecMag(MC) = mag
      WRITE(iFile,*)MC, mag

      Call ClassificationMatrix(L, S, C, PBC, degen, vecEps)
      MC = 2
      DO WHILE (MC.LE.MCS)
        spinType = ClassChoice(vecTransProb, degen)
        g_l = degen(spinType)           
        spinK = INT(g_l*RAND()) + 1

        Call CalculateDeltat(vecTransProb , degen, Delta_t, Q)

        iMax = MIN(MC + Delta_t, MCS)
        DO i = MC + 1, iMax
          WRITE(iFile,*)i, DABS(mag)
          vecMag(i) = DABS(mag)
        ENDDO
 
        IF (DABS(Q - 1.d0).LT.10.d0**(-13.d0)) THEN
            DO i = MC + 1, MCS
              vecMag(i) = DABS(mag)
              WRITE(iFile,*)i, DABS(mag)
            ENDDO
          EXIT
          ENDIF

        Call FlipChangeClass(N, L, spinK, degen, spinType, C, 
     &    S, vecEps, mag, PBC)

      MC = MC + Delta_t
      ENDDO
 
      WRITE(iFile,*)
      WRITE(iFile,*)
      RETURN
      END

C-----Computes de time a spin will take before it turns-----------------
      SUBROUTINE CalculateDeltat(vecTransProb, degen, Delta_t, Q)
      !-----Variables---------------------------------------------------
      !INPUTS
      !vecTransProb -> Transition probabilities for each type
      !degen -> Vector of degeneration for each type of transition
      !OUTPUTS
      !Delta_t -> Time spin will take before it turns
      !Q -> Number of arrows occupying each one of the different possible states
      !-----------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER i, degen(5), Delta_t
      DOUBLE PRECISION vecTransProb(5), suma, Q

      suma = 0.d0
      DO i = 1, 5
        suma = suma + vecTransProb(i)*degen(i)
      ENDDO
      Q = 1.d0 - suma

      Delta_t = 1 + NINT(DLOG(DBLE(RAND()))/(DLOG(Q)))

      RETURN
      END

C-----Gets the index of the N disks randomly----------------------------
      SUBROUTINE FlipChangeClass(N, L, spinK, degen, spinType,
     &  C, S, vecEps, mag, PBC) 
      !-----Variables---------------------------------------------------
      !index_mov---> Npart vector with randomly selected indices
      !1 --> 2; 2 --> 1; 3 --> 4; 4 --> 3; 5 --> 5
      !-----------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER i, j, L, degen(5), spinType, spinK,
     &  C(1:L,1:L), vecEps(5), energy, S(1:L, 1:L), N, k, 
     &  PBC(0:L+1)
      DOUBLE PRECISION mag, pi
      pi = 4*DATAN(1.D0)

      k = 0
      i = 1
      j = 1
      DO WHILE (k.LT.spinK)
        energy = C(i, j)

        IF (energy.EQ.(vecEps(spinType))) THEN
          k = k + 1
          IF (k.LT.spinK) THEN
            IF (j.LT.L) THEN
              j = j + 1
            ELSE
              i = i + 1
              j = 1
            ENDIF
           ENDIF
        ELSE
          IF (j.LT.L) THEN
            j = j + 1
          ELSE
            i = i + 1
            j = 1
          ENDIF
        ENDIF

      ENDDO

      S(i,j) = -S(i,j)
      mag = mag + (2*S(i,j))/DBLE(N)

      CALL ClassificationMatrix(L, S, C, PBC, degen, vecEps)
        
      RETURN
      END