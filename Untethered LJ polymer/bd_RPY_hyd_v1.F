#include "fintrf.h"
C D = bd_RPY_hyd_v1(np,n3,a,L,mu,x,y,z);
C np = number of beads, n3 = 3*np, a= bead radius, L= bead diameter
C mu = viscosity
C x,y,z = positon
C D = hydrodynamic tensor
C Compilied in matlab "mex bd_RPY_hyd_v1.F"  
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
      mwPointer a,L,mu,x,y,z
      mwPointer D

C     Array Information
      mwPointer m,n,o
      
C-----------------------------------------------------------------------
C     Get the size of the input array.
      
      m = mxGetM(prhs(6))
      n = mxGetN(prhs(6))
      o = INT(3*m)

C     Create Fortran array from the input argument

      a  = mxGetPr(prhs(3))
      L  = mxGetPr(prhs(4))
      mu = mxGetPr(prhs(5))
      x  = mxGetPr(prhs(6))
      y  = mxGetPr(prhs(7))
      z  = mxGetPr(prhs(8))
      
      plhs(1) = mxCreateDoubleMatrix(o,o,0)
      D    = mxGetPr(plhs(1))



C     Call the computational subroutine.
      call for(%VAL(D),%VAL(a),%VAL(L),%VAL(mu),
     +              %VAL(x),%VAL(y),%VAL(z),m,o)
      RETURN
      END
C-----------------------------------------------------------------------
C     Computational routine

      SUBROUTINE for(D,a,L,mu,x,y,z,m,o)
      implicit none

      mwSize m,o,I,J,IC,JC
      real*8 RD(o,o)
      real*8 D(o,o), x(m), y(m), z(m),a,L,mu
      real*8 RXIJ, RYIJ, RZIJ, RXIJSQ, RYIJSQ, RZIJSQ
      real*8 RIJSQ, RIJ, R6I, R12I, RXI, RYI,RZI
      real*8 C0, C1, C2, C3, C4
      real*8 FXBF, FYBF, FZBF, DIF, FIJ, RCUT, RMAX, A2, A6
      real*8 IZIJ, IZIJSQ
      real*8 A4XX, A4YY, A5ZZ, A6XY, A6YX, A7XZ, A7YZ, A8ZX, A8ZY, A7TM 
      real*8 R2, R3, R5, R7, AR5, AR7, ZZIJ, AAA 
      real*8 AZI, L6, CONST, CONST2
      real*8 IIJ, IIJSQ, PI

      A2   = a**2
      L6   = L**6
      D    = 0.0
      PI   = 3.141592653
     
      RD   = 0
     
      do I=1,m-1
        
    	IC = 3*(I-1)+1
        RXI =  x(I)
      	RYI =  y(I)
      	RZI =  z(I)
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
      	
        if (RIJSQ .gt. 4*a) then 
            C0 = 3 *  a/(4*RIJ)
            C1 = 1 + 2*A2/(3*RIJSQ)
            C2 = 1 - 2*A2/RIJSQ
            C3 = C0 * C1
            C4 = C0 * C2
            RD( IC  , JC   ) = C3 + C4 * RXIJSQ/RIJSQ
            RD( IC+1, JC+1 ) = C3 + C4 * RYIJSQ/RIJSQ
            RD( IC+2, JC+2 ) = C3 + C4 * RZIJSQ/RIJSQ
            RD( IC  , JC+1 ) = C4 * RXIJ * RYIJ/RIJSQ
            RD( IC  , JC+2 ) = C4 * RXIJ * RZIJ/RIJSQ
            RD( IC+1, JC+2 ) = C4 * RYIJ * RZIJ/RIJSQ
            RD( IC+1, JC   ) = RD( IC  , JC+1 )
            RD( IC+2, JC   ) = RD( IC  , JC+2 )
            RD( IC+2, JC+1 ) = RD( IC+1, JC+2 )  
        ELSE
            CONST=1.0-(9.0*RIJ/(a*32.0))
            CONST2=3.0/(32.0*a*RIJ)
            RD( IC  , JC   ) = CONST + CONST2 * RXIJSQ
            RD( IC+1, JC+1 ) = CONST + CONST2 * RYIJSQ
            RD( IC+2, JC+2 ) = CONST + CONST2 * RZIJSQ
            RD( IC  , JC+1 ) = CONST2 * RXIJ * RYIJ
            RD( IC  , JC+2 ) = CONST2 * RXIJ * RZIJ
            RD( IC+1, JC+2 ) = CONST2 * RYIJ * RZIJ
            RD( IC+1, JC   ) = RD( IC  , JC+1 )
            RD( IC+2, JC   ) = RD( IC  , JC+2 )
            RD( IC+2, JC+1 ) = RD( IC+1, JC+2 )  
        endif 
      end do 

      end do

      RD = (RD + TRANSPOSE(RD))
      D(:,:) = RD(:,:)

C     On diagonal componets of blake tensor and 
C     Force Interaction with surface
      do I = 1, m
           IC = 3 * ( I - 1 ) + 1
           D( IC  , IC   ) = 1
           D( IC+1, IC+1 ) = 1
           D( IC+2, IC+2 ) = 1
      end do
      RETURN
      END