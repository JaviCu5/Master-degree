C2345678901234567890123456789012345678901234567890123456789012
      PROGRAM practice_1
      IMPLICIT NONE
      INTEGER ChemSp, R, i, nbox, nwrite, write, N
      PARAMETER (ChemSp=2, R=4, nbox=40, N=1000)
      DOUBLE PRECISION params(R), fm, fp, Specie(ChemSp),
     &  Vol, t_fin,toConc,h, mstat(N), pstat(N), xmRNA(nbox),
     & mprob(nbox), merr(nbox), xprot(nbox), pprob(nbox),
     & perr(nbox),mMean,mdesv,pMean,pdesv,Mean,F1,FnegDet,
     & FnegGill
      EXTERNAL F1,FnegDet,FnegGill

      COMMON /dummyVar/ t_fin, toConc, nwrite, write, Vol,h

      !Some parameters
      fm = 1.0d0
      fp = 1.0d0
C The params vector contains alpham, deltam, alphap, deltap in this order
      params(1) = 100.0d0*fm
      params(2) = 1.0d0*fm
      params(3) = 10.0d0*fp
      params(4) = 0.1d0*fp

      Vol = 1d0 !um^3  this is the volume used on Gillespie algorithm
      toConc = 1/(Vol*0.6022) !the parameter for convertin molecules to concentration
      nwrite = 100 !write 1 date on 100
      h=0.001d0 !the integration step in RK4 subroutine

      OPEN(10, file="Gillespie.dat")
      OPEN(11, file="Deterministic.dat")
      OPEN(12, file="Histograms.dat")

C-----------------------Default case--------------------------
      t_fin = 50.d0 !final time for the deterministic and stochastic simulation in the default case

      write=1
      WRITE(11,*)"# INDEX 0 time evolution IC 0, 0"
C     Deterministic subroutine 
      CALL DeterministicSimulation(ChemSp, R, params,F1)
      DO i = 1, 2
        CALL srand(10*i)
        !Initial conditions
        Specie(1)=0.0d0
        Specie(2)=0.0d0
        WRITE(10,*)"INDEX ", i-1, " time evolution IC 0, 0"
C       Stochastic subroutine
        CALL GillespieSimulation(params,Specie(1),Specie(2),R,
     &  F1)
      ENDDO

      t_fin = 10.d0 !final time for the confections of histograms
      write = 0 !parameter for not writting the data
      DO i = 1, N
        CALL srand(10*i)
        !Initial conditions, deterministic steady values
        Specie(1)=100.0d0/toConc
        Specie(2)=9925.0d0/toConc
        CALL GillespieSimulation(params,Specie(1),Specie(2),R,
     &  F1)
C       saving the steady state for mRNA and protein concentrations
        mstat(i) = Specie(1)*toConc
        pstat(i) = Specie(2)*toConc
      ENDDO
C     mRNA concentration steady state histogram
      CALL histograma(N,mstat,60.d0,150.d0,nbox,xmRNA,
     & mprob, merr) !60,140  92,108
C     protein concentration steady state histogram
      CALL histograma(N,pstat,9000.d0,11200.d0,nbox,xprot,
     & pprob, perr)

C     Compute the gausian parameters for comparison with theory
      mMean = Mean(mstat, N)
      mdesv = sqrt(mMean*toConc)
      pMean = Mean(pstat, N)
      pdesv = sqrt(pMean*toConc)

      WRITE(*,*)"Mean and desv: m, pr",mMean,mdesv,pMean,pdesv
      WRITE(10,*)"INDEX 0 prob density mRNA and prot"
      DO i = 1, nbox
       WRITE(12,*)xmRNA(i),mprob(i), merr(i),xprot(i),
     & pprob(i), perr(i)
      ENDDO
      WRITE(12,*)
      WRITE(12,*)


C-------Effect of transcription and translation rates---------
      fm = 1.0d0
      fp = 1.0d0
C The params vector contains alpham, deltam, alphap, deltap in this order
      params(1) = 100.d0*10.d0*fm
      params(2) = 1.0d0*fm
      params(3) = 10.0d0/10.d0*fp
      params(4) = 0.1d0*fp

      write=1
      t_fin=70.d0
      WRITE(11,*)"# INDEX 1 time evolution rates effect"
      CALL DeterministicSimulation(ChemSp, R, params,F1)
      DO i = 1, 2
        CALL srand(10*i)
        !Initial conditions
        Specie(1)=0.0d0
        Specie(2)=0.0d0
        WRITE(10,*)"INDEX ",i+1, " time evol rates effect"
        CALL GillespieSimulation(params,Specie(1),Specie(2),R,
     &  F1)
      ENDDO

      t_fin = 10.d0
      write = 0
      DO i = 1, N
        CALL srand(10*i)
        !Initial conditions
        Specie(1)=1000.0d0/toConc
        Specie(2)=9989.868d0/toConc
        CALL GillespieSimulation(params,Specie(1),Specie(2),R,
     &  F1)
        mstat(i) = Specie(1)*toConc
        pstat(i) = Specie(2)*toConc
      ENDDO

      CALL histograma(N,mstat,800.d0,1200.d0,nbox,xmRNA,
     & mprob, merr) 

      CALL histograma(N,pstat,9200.d0,10800.d0,nbox,xprot,
     & pprob, perr)

      mMean = Mean(mstat, N)
      mdesv = sqrt(mMean*toConc)
      pMean = Mean(pstat, N)
      pdesv = sqrt(pMean*toConc)

      WRITE(*,*)"Mean and desv rates: m, pr",mMean,mdesv,
     &  pMean,pdesv
      WRITE(10,*)"INDEX 1 prob density mRNA and prot rates"
      DO i = 1, nbox
       WRITE(12,*)xmRNA(i),mprob(i), merr(i),xprot(i),
     & pprob(i), perr(i)
      ENDDO
      WRITE(12,*)
      WRITE(12,*)

C----------------Effect of negative feedback------------------
      fm = 1.0d0
      fp = 1.0d0
C The params vector contains alpham, deltam, alphap, deltap in this order
      params(1) = 10100.d0*fm
      params(2) = 1.0d0*fm
      params(3) = 10.0d0*fp
      params(4) = 0.1d0*fp

      write=1
      t_fin=30.d0
      WRITE(11,*)"# INDEX 2 time evolution negative feedback"
      CALL DeterministicSimulation(ChemSp, R, params,FnegDet)
      DO i = 1, 2
        CALL srand(10*i)
        !Initial conditions
        Specie(1)=0.0d0
        Specie(2)=0.0d0
        WRITE(10,*)"INDEX ",i+3, " time evol negative feedback"
        CALL GillespieSimulation(params,Specie(1),Specie(2),R,
     &  FnegGill)
      ENDDO

      t_fin = 10.d0
      write = 0
      DO i = 1, N
        CALL srand(10*i)
        !Initial conditions
        Specie(1)=458.544d0/toConc
        Specie(2)=4585.434d0/toConc
        CALL GillespieSimulation(params,Specie(1),Specie(2),R,
     &  FnegGill)
        mstat(i) = Specie(1)*toConc
        pstat(i) = Specie(2)*toConc
      ENDDO

      CALL histograma(N,mstat,60.d0,150.d0,nbox,xmRNA,
     & mprob, merr) !60,140  92,108

      CALL histograma(N,pstat,9000.d0,11200.d0,nbox,xprot,
     & pprob, perr)

      mMean = Mean(mstat, N)
      mdesv = sqrt(mMean*toConc)
      pMean = Mean(pstat, N)
      pdesv = sqrt(pMean*toConc)

      WRITE(*,*)"Mean and desv rates: m, pr",mMean,mdesv,
     &  pMean,pdesv
      WRITE(10,*)"INDEX 1 prob density mRNA and prot rates"
      DO i = 1, nbox
       WRITE(12,*)xmRNA(i),mprob(i), merr(i),xprot(i),
     & pprob(i), perr(i)
      ENDDO
      WRITE(12,*)
      WRITE(12,*)
      
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      
      END

C-------------------------------------------------------------
C-------------------------------------------------------------
C--------------------------FUNCTIONS--------------------------
C-------------------------------------------------------------
C-------------------------------------------------------------

      INTEGER FUNCTION ComputeNu(r2, a_mu, R)
      IMPLICIT NONE
      INTEGER R
      DOUBLE PRECISION a_mu(R), aux, r2
      !This function computes the nu on the Gillespie algorithm in order 
      !to know which reaction is activated in each MC step
      ComputeNu = 1
      aux = a_mu(ComputeNu)

      DO WHILE (aux.lt.r2)
            ComputeNu = ComputeNu + 1
            aux = aux + a_mu(ComputeNu)
      ENDDO
      END

C-------------------------------------------------------------
      DOUBLE PRECISION FUNCTION Mean(vect, ndat)
      IMPLICIT NONE
      INTEGER ndat, i
      DOUBLE PRECISION vect(ndat)
      !This function computes the mean of a vvector of values

      Mean = 0.d0
      DO i = 1, ndat
            Mean = Mean + Vect(i)
      ENDDO

      Mean = Mean/ndat
      END

C-------------------------------------------------------------
      DOUBLE PRECISION FUNCTION F1(p)
      IMPLICIT NONE 
      DOUBLE PRECISION p
      !This is how affect the protein concentration to the mRNA generation
      F1 = 1
      END FUNCTION F1

C-------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FnegDet(p)
      IMPLICIT NONE 
      DOUBLE PRECISION p,k
      !This is other function for the influence of the protein concentration to the 
      !mRNA generation, special for the molecular Gillespie algorithm
      k=1000
      FnegDet = (k*k)/(k*k+p*p)
      END FUNCTION FnegDet

C-------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FnegGill(p)
      IMPLICIT NONE 
      DOUBLE PRECISION p,K
      !This is other function for the influence of the protein concentration to the 
      !mRNA generation, special for deterministic RungeKutta4.
      K=602
      FnegGill = (K*K)/(K*K+p*p)
      END FUNCTION FnegGill

C-------------------------------------------------------------
C-------------------------------------------------------------
C--------------------------SUBROUTINES------------------------
C-------------------------------------------------------------
C-------------------------------------------------------------

      SUBROUTINE Compute_a_mu(params,M,P,a_mu,a_0,R,Vol,func)
      IMPLICIT NONE
      INTEGER R, i
      DOUBLE PRECISION a_mu(R),a_0,params(R),M,P,Vol,func
      EXTERNAL func
      !This computes the propensities for given protein influence
      !It also computes the accumulative propensity 
      
      a_mu(1) = Vol*params(1)*0.6022d0*func(P)
      a_mu(2) = params(2)*M
      a_mu(3) = params(3)*M
      a_mu(4) = params(4)*P

      a_0 = 0d0
      DO i = 1, R
            a_0 = a_0 + a_mu(i)
      ENDDO

      END

C-------------------------------------------------------------
      SUBROUTINE UpdateSpecies(M, P, nu)
      IMPLICIT NONE
      INTEGER nu
      DOUBLE PRECISION M, P
      !This adds/substracs a molecule of mRNA/Protein for each MC step in the
      !Gilispie algorithm for a given nu (reaction)
      IF (nu.eq.1) THEN
            M = M + 1
      ELSE IF (nu.eq.2) THEN
            M = M - 1
      ELSE IF (nu.eq.3) THEN
            P = P + 1
      ELSE
            P = P - 1
      ENDIF
      END

C-------------------------------------------------------------
      SUBROUTINE MCStep(a_mu, r2, M, P, R)
      !This Monte Carlo step choose whitch reaction activates and update the chemical specie related
      IMPLICIT NONE
      INTEGER R, nu, ComputeNu
      DOUBLE PRECISION a_mu(R), r2, M, P

      nu = ComputeNu(r2, a_mu, R) !Fins the reaction activated
      CALL UpdateSpecies(M, P, nu) !Updates the count molecules

      END

C-------------------------------------------------------------
      SUBROUTINE GillespieSimulation(params, M, P, R,func)
      !It computes a stochastic Gillespie simulation from initial conditions
      !to a given final time (t_fin)
      IMPLICIT NONE
      INTEGER R, write, j, nwrite
      DOUBLE PRECISION a_mu(R), a_0, tau, toConc, h,func,
     &  r2, M, P, t, t_fin, params(R), Vol
      EXTERNAL func

      COMMON /dummyVar/ t_fin, toConc, nwrite, write, Vol,h
      
      t = 0.d0
      j=0
      DO WHILE (t.le.t_fin)
            IF(write.eq.1)THEN
            IF(mod(j, nwrite).eq.0)then
                  WRITE(10,*)t, M*toConc, P*toConc
            ENDIF
            ENDIF
            CALL Compute_a_mu(params,M,P,a_mu,a_0,R,Vol,
     &       func) !computes the propensities
     
            tau = - log(1-RAND())/a_0 !Computes the expected time for the actual MC step
            r2 = a_0*RAND() !Computes the random accumulated propensities 
            CALL MCStep(a_mu, r2, M, P, R) !Calls the MC step to update the system

            t = t + tau !Computes the actal time after the MC step
            j=j+1 !Just for knowing which data we save
      ENDDO

      IF(write.eq.1)THEN !For the gnuplot to separate graphisc
            WRITE(10,*)
            WRITE(10,*)
      ENDIF

      END

C-------------------------------------------------------------
      SUBROUTINE DeterministicSimulation(nequ, R, params,func)
      !It computes a deterministic simulation from initial conditions
      !to a given final time (t_fin)
      IMPLICIT NONE
      INTEGER npas, nwrite, nequ, i, j, R, write
      DOUBLE PRECISION t_fin, h, t, y(nequ), yout(nequ),
     &  params(R), toConc, Vol, func
      EXTERNAL func

      COMMON /dummyVar/ t_fin, toConc, nwrite, write, Vol,h

      npas = int(t_fin/h) !How many RK4 steps we will need

      do j=1,npas
           
      call RK4(func,nequ,t,h,y,yout,R,params)
ccccc    update independent variable
      t=dble(j)*h

ccccc    write data 
      IF((write.eq.1).and.(mod(j,nwrite).eq.0))then
        WRITE(11,*)t,yout(1),yout(2)
      ENDIF
     
ccccc    update variables for next step
      DO i=1,nequ
        y(i)=yout(i)
      ENDDO
      ENDDO
      WRITE(11,*)
      WRITE(11,*)
      END

C-------------------------------------------------------------
      SUBROUTINE histograma(ndat,xdata,xa,xb,nbox,xhisto,
     & histo,errhisto)
      IMPLICIT NONE
      INTEGER i, k, nbox, ndat, ntot
      DOUBLE PRECISION xa, xb, xdata(ndat), xhisto(nbox),
     & histo(nbox), errhisto(nbox), deltax
       
      deltax = (xb - xa)/dble(nbox)
      DO k = 1, nbox 
            xhisto(k)=xa-0.5d0*deltax+dble(k)*deltax
            histo(k)=0.0d0
      ENDDO
      ntot=0
      
      DO i = 1, ndat
      IF((xdata(i).gt.xb).or.(xdata(i).lt.xa))THEN
            WRITE(6,*)'data out of xa-xb range:'
            WRITE(6,*)xdata(i)
      ELSE
            k=int((xdata(i)-xa)/deltax)+1  
            histo(k)=histo(k)+1.0d0
            ntot=ntot+1
      ENDIF
      ENDDO
      WRITE(6,*)'the histogram only includes',ntot,'data from the total'
                 
c constructs probability density and error (assuming binomial)
      DO k=1,nbox
            histo(k)=histo(k)/(dble(ntot))
            errhisto(k)=dsqrt(histo(k)*(1.0d0-histo(k))/dble(ntot))/
     &        deltax
            histo(k)=histo(k)/deltax
      ENDDO
       
      RETURN
      END

C-------------------------------------------------------------
      SUBROUTINE RK4(func,nequ, x, h, y, yout, R, params)
      IMPLICIT NONE 
      INTEGER nequ, R
      DOUBLE PRECISION y(nequ), yout(nequ), dy1(nequ), 
     &  dy2(nequ), dy3(nequ), dy4(nequ), x1, x3,
     &  y1(nequ), y2(nequ), y3(nequ), x, h, params(R),func
      EXTERNAL func
      !A typical order 4vRungeKutta soubroutine for a given system
      
      CALL DFUNCTIONS(func,NEQU, y, dy1, R, params)
      x1=x+h*0.5d0
      y1=y+h*0.5d0*dy1  
      CALL DFUNCTIONS(func,NEQU, y1, dy2, R, params)
      y2=y+h*0.5d0*dy2  
      CALL DFUNCTIONS(func,NEQU, y2, dy3, R, params)
      x3=x+h
      y3=y+h*dy3
      CALL DFUNCTIONS(func,NEQU, y3, dy4, R, params)
      
      yout=y+h*(dy1+2.0d0*(dy2+dy3)+dy4)/6.0d0
      RETURN
      END

C-------------------------------------------------------------
      SUBROUTINE DFUNCTIONS(func,NEQU, y, dy, R, params)
      IMPLICIT NONE
      INTEGER NEQU, R
      DOUBLE PRECISION y(NEQU),dy(NEQU),params(R),func
      !Descriptive ordinary differential equations for the system

        dy(1)=params(1)*func(y(2))-params(2)*y(1)
        dy(2)=params(3)*y(1)-params(4)*y(2)

      RETURN
      END
