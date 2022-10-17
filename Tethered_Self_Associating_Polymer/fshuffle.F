#include "fintrf.h"
C-----------------------------------------------------------------------
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
 
C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs,intbond

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetM, mxGetN

C     Pointers to Input and Output mxArrays:
      mwPointer In
      mwPointer R

C     Array Information
      mwPointer m,n
      
C-----------------------------------------------------------------------
C     Get the size of the input array.
      
C     Create Fortran array from the input argument

      In  = mxGetPr(prhs(1))

      plhs(1) = mxCreateDoubleMatrix(1,1,0)
      R = mxGetPr(plhs(1))
   
C     Call the computational subroutine.
      call shuffle(%VAL(In),%VAL(R))
      RETURN
      END
C-----------------------------------------------------------------------
C     Computational routine

      SUBROUTINE shuffle(In,R)
      implicit none
      Integer seed,I
      INTEGER, ALLOCATABLE :: new (:), old(:)
      REAL*8 R, In
      seed=int(In)
      CALL RANDOM_SEED ( )  
                      
      CALL RANDOM_SEED (SIZE = I)  
                             
      ALLOCATE (new(I))
      ALLOCATE (old(I))
      CALL RANDOM_SEED (GET=old(1:I)) 
      new = seed
      CALL RANDOM_SEED (PUT=new(1:I)) 
      CALL Random_Number(R)                	     
      RETURN    
      END