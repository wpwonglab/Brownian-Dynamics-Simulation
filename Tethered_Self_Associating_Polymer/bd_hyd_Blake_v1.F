#include "fintrf.h"
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
      real*8 RD(o,o), ID(o,o), BU(o,o), BL(o,o)
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
      ID   = 0
      BU   = 0
      BL   = 0
     
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

            IZIJ   = RZI + z(J)
            IZIJSQ = IZIJ**2
       	    IIJSQ  = RXIJSQ + RYIJSQ + IZIJSQ
            IIJ    = SQRT(IIJSQ)

            C0 = 3 * a/(4*IIJ)  
            C1 = 1 + 2*A2/(3*IIJSQ)
            C2 = 1 - 2*A2/IIJSQ
            C3 = C0 * C1
            C4 = C0 * C2
            ID( IC  , JC   ) = C3 + C4 * RXIJSQ/IIJSQ
            ID( IC+1, JC+1 ) = C3 + C4 * RYIJSQ/IIJSQ
            ID( IC+2, JC+2 ) = C3 + C4 * IZIJSQ/IIJSQ
            ID( IC  , JC+1 ) = C4 * RXIJ * RYIJ/IIJSQ
            ID( IC  , JC+2 ) = C4 * RXIJ * IZIJ/IIJSQ
            ID( IC+1, JC+2 ) = C4 * RYIJ * IZIJ/IIJSQ
            ID( IC+1, JC   ) = ID( IC  , JC+1 )
            ID( IC+2, JC   ) = ID( IC  , JC+2 )
            ID( IC+2, JC+1 ) = ID( IC+1, JC+2 )

            R2   = IIJ**2
            R3   = R2*IIJ
            R5   = R3*R2
            R7   = R5*R2
            AR5  = a**2/R5
            AR7  = AR5/R2
            AAA  = RZI*z(J)/R3

            A4XX = AAA*(3*RXIJSQ/R2-1)+AR5*IZIJSQ*(1-5*RXIJSQ/R2)
            A4YY = AAA*(3*RYIJSQ/R2-1)+AR5*IZIJSQ*(1-5*RYIJSQ/R2)
            A5ZZ = AAA*(1-3*IZIJSQ/R2)-AR5*IZIJSQ*(3-5*IZIJSQ/R2)
            A6XY = 3*RZI*z(J)*RXIJ*RYIJ/R5-5*AR7*RXIJ*RYIJ*IZIJSQ
            A6YX = A6XY
			
            A7TM = 2-5*IZIJSQ/R2
            A7XZ = z(J)*RXIJ/R3*(1-3*z(I)*IZIJ/R2)-AR5*RXIJ*IZIJ*A7TM
            A7YZ = z(J)*RYIJ/R3*(1-3*z(I)*IZIJ/R2)-AR5*RYIJ*IZIJ*A7TM
            A8ZX = z(J)*RXIJ/R3*(1+3*z(I)*IZIJ/R2)-5*AR7*RXIJ*IZIJ**3
            A8ZY = z(J)*RYIJ/R3*(1+3*z(I)*IZIJ/R2)-5*AR7*RYIJ*IZIJ**3

            BU( IC  , JC   ) = A4XX
            BU( IC+1, JC+1 ) = A4YY
            BU( IC+2, JC+2 ) = A5ZZ
            BU( IC  , JC+1 ) = A6XY
            BU( IC  , JC+2 ) = A7XZ
            BU( IC+1, JC+2 ) = A7YZ
            BU( IC+1, JC   ) = A6YX
            BU( IC+2, JC   ) = A8ZX
            BU( IC+2, JC+1 ) = A8ZY



      end do 

      end do

      RD = (RD + TRANSPOSE(RD))
      ID = (ID + TRANSPOSE(ID))
      BU = (BU + TRANSPOSE(BU))
      BU = 3*a/2*BU
        
      D(:,:) = RD(:,:)-ID(:,:)+BU(:,:)

C     On diagonal componets of blake tensor and 
C     Force Interaction with surface
      do I = 1, m
           IC = 3 * ( I - 1 ) + 1
           AZI = a/z(I)
           D( IC  , IC   ) = (1-0.5625*AZI+0.125*AZI**3)
           D( IC+1, IC+1 ) = D(IC, IC)
           D( IC+2, IC+2 ) = (1-1.1250*AZI+0.500*AZI**3)
      end do
      RETURN
      END