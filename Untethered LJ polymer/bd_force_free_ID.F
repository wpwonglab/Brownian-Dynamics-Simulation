#include "fintrf.h"
C Calculate forces for polymer 
C [~, TFX, TFY, TFZ, ~] =bd_force_free_ID(np,n3,u,a,K,L,x,y,z,mu) 
C np = number of beads, n3 = 3*np, u =lennard jones potential
C a = length of bead, K = spring constant, L = bead diameter
C xyz= position, mu =  viscosity
C TFX, TFY,TFZ -> forces 
C Compilied in matlab "mex bd_force_free_ID.F"   
C-----------------------------------------------------------------------
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
 
C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetM, mxGetN

C     Pointers to Input and Output mxArrays:
      mwPointer U,a,K,L,x,y,z,mu
      mwPointer ET,TFX,TFY,TFZ

C     Array Information
      mwPointer m,n,o
      
C-----------------------------------------------------------------------
C     Get the size of the input array.
      
      m = mxGetM(prhs(7))
      n = mxGetN(prhs(7))
      o = INT(3*m)

C     Create Fortran array from the input argument

      u  = mxGetPr(prhs(3))
      a  = mxGetPr(prhs(4))
      K  = mxGetPr(prhs(5))
      L  = mxGetPr(prhs(6))
      x  = mxGetPr(prhs(7))
      y  = mxGetPr(prhs(8))
      z  = mxGetPr(prhs(9))
      mu = mxGetPr(prhs(10))



      plhs(1) = mxCreateDoubleMatrix(1,3,0)
      plhs(2) = mxCreateDoubleMatrix(m,n,0)
      plhs(3) = mxCreateDoubleMatrix(m,n,0)
      plhs(4) = mxCreateDoubleMatrix(m,n,0)
      plhs(5) = mxCreateDoubleMatrix(m,n,0)
      plhs(6) = mxCreateDoubleMatrix(m,m,0)
      plhs(7) = mxCreateDoubleMatrix(m,m,0)




      ET   = mxGetPr(plhs(1))
      TFX  = mxGetPr(plhs(2))
      TFY  = mxGetPr(plhs(3))
      TFZ  = mxGetPr(plhs(4))



C     Call the computational subroutine.
      call for(%VAL(ET),%VAL(TFX),%VAL(TFY),%VAL(TFZ),
     +              %VAL(U),%VAL(a),%VAL(K),%VAL(L),%VAL(mu),
     +              %VAL(x),%VAL(y),%VAL(z),m,o,)
      RETURN
      END
C-----------------------------------------------------------------------
C     Computational routine

      SUBROUTINE for(ET,TFX,TFY,TFZ,u,a,K,L,mu,x,y,z,m,o)
      implicit none

      mwSize m,o,I,J,IC,JC
      real*8 ET(1,3), TFX(m), TFY(m), TFZ(m)
      real*8 u, a, K, L, mu, x(m), y(m), z(m), A2, RMAX
      real*8 SURF(m)
      real*8 ELJ, EBF,ESURF
      real*8 FX(m), FY(m), FZ(m), BX(m), BY(m), BZ(m)
      real*8 RXI, RYI, RZI, FXI, FYI, FZI, BXI, BYI, BZI
      real*8 RXIJ, RYIJ, RZIJ, RXIJSQ, RYIJSQ, RZIJSQ
      real*8 RIJSQ, RIJ, R6I, R12I, FENE
      real*8 FXIJ, FYIJ, FZIJ, L6
      real*8 FXBF, FYBF, FZBF, DIF, FIJ, RCUT
      real*8 IIJ, IIJSQ, PI
      real*8 UTX,UTY,UTZ,TRAPX,TRAPY,TRAPZ,ETRAP, FDIF

      RCUT = 12*L
      A2   = a**2
      L6   = L**6
      PI   = 3.141592653
      RMAX =  1.5 * 2 * a /2**(1/6)
     

      FX  = 0.0
      FY  = 0.0
      FZ  = 0.0
      BX  = 0.0
      BY  = 0.0
      BZ  = 0.0 

      TFX = 0.0
      TFY = 0.0
      TFZ = 0.0
      ELJ = 0.0
      EBF = 0.0	 

        
      do I=1,m-1

    	IC = 3*(I-1)+1
    
      	RXI =  x(I)
      	RYI =  y(I)
      	RZI =  z(I)
      	FXI = FX(I)
      	FYI = FY(I)
      	FZI = FZ(I)
      	BXI = BX(I)
      	BYI = BY(I)
    	BZI = BZ(I)
    
    	do J=I+1,m

            JC = 3*(J-1)+1

            RXIJ = RXI - x(J)
            RYIJ = RYI - y(J)
            RZIJ = RZI - z(J)
            RXIJSQ = RXIJ**2
            RYIJSQ = RYIJ**2
            RZIJSQ = RZIJ**2
       	    RIJSQ = RXIJSQ + RYIJSQ + RZIJSQ
            RIJ   = SQRT(RIJSQ)

C     Calculate LJ forces
      if (RIJ .LT. RCUT ) then
            R6I   = (L/RIJ)**6
            R12I  = R6I**2
            ELJ   = ELJ + (R12I - 2*R6I)
            FIJ   = 12*(R12I-R6I)/RIJSQ
            FXIJ  = RXIJ * FIJ
            FYIJ  = RYIJ * FIJ
            FZIJ  = RZIJ * FIJ
            FXI   = FXI + FXIJ
            FYI   = FYI + FYIJ
            FZI   = FZI + FZIJ
            FX(J) = FX(J) - FXIJ
            FY(J) = FY(J) - FYIJ
            FZ(J) = FZ(J) - FZIJ
      endif
      



C     Calculate Bond forces
      if(J .eq. (I+1)) then

            FDIF  = RIJ-L 
            FENE  = 1 - (FDIF/RMAX)**2
            DIF   = FDIF/FENE
            FXBF  = DIF * RXIJ / RIJ
            FYBF  = DIF * RYIJ / RIJ
            FZBF  = DIF * RZIJ / RIJ
            BXI   = BXI - FXBF
            BYI   = BYI - FYBF
            BZI   = BZI - FZBF
            BX(J) = FXBF
            BY(J) = FYBF
            BZ(J) = FZBF
            EBF = EBF + LOG(FENE)
       endif


      end do 
C     LJ force   
 
      FX(I) = FXI
      FY(I) = FYI
      FZ(I) = FZI

C    Spring force for adjacent bonds

      BX(I) = BXI
      BY(I) = BYI
      BZ(I) = BZI  
   

      end do


C       Sum Total energies and Forces 
        ET(1,1) = u*ELJ 
        ET(1,2) = K/2*(EBF)
        ET(1,3) = ET(1,1)+ET(1,2)
        TFX = u*FX + K*BX
        TFY = u*FY + K*BY
        TFZ = u*FZ + K*BZ




      RETURN
      END