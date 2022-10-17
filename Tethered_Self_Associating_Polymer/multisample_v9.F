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
      mwPointer np,a,vollap,fbreak,rmx,delt,kon,bcut,bond
      mwPointer rmxout

C     Array Information
      mwPointer m,n
      
C-----------------------------------------------------------------------
C     Get the size of the input array.
      
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
C     Create Fortran array from the input argument

      np  = mxGetPr(prhs(1))
      a  = mxGetPr(prhs(2))
      vollap  = mxGetPr(prhs(3))
      fbreak  = mxGetPr(prhs(4))
      rmx  = mxGetPr(prhs(5))
      delt = mxGetPr(prhs(6))
      kon  = mxGetPr(prhs(7))
      bcut = mxGetPr(prhs(8))
      bond = mxGetPr(prhs(9))
        

      plhs(1) = mxCreateDoubleMatrix(m,m,0)
        
      rmxout = mxGetPr(plhs(1))


C     Call the computational subroutine.
      call sample(%VAL(np),%VAL(a),%VAL(vollap),%VAL(fbreak),%VAL(rmx),
     +              %VAL(delt),%VAL(rmxout),m,n,
     +              %VAL(kon),%VAL(bcut),%VAL(bond))
      RETURN
      END
C-----------------------------------------------------------------------
C     Computational routine

      SUBROUTINE sample(np,a,vollap,fbreak,rmx,delt,
     +             rmxout,m,n,kon,bcut,bond)
      implicit none

      INTEGER m,n, I, J, K, C1, C2, B1, B2, r1, r2, A1, A2, A3, A4
      INTEGER C3,C4, bint, MS(m,10)
      real*8 a,Eb, Eo, delt,isbonded, r3, r4,np, Eub, vhem,PI,vin
      real*8 vollap(m,m), fbreak(m,m), rmx(m,m),rmxout(m,m)
      real*8 array_1(m), array_2(m),prob(m,m), probbreak(m,m)
      real*8 koff,xt, koffav,RMAX,bcut,kon,bond
      real*8 seed
      rmxout = rmx
      bint = int (bond)

      PI   = 3.141592653
      vhem = (4.0 * PI *((bcut/2)*a)**3)/3.0
      vin = (4.0 * PI *(a*0.94)**3)/3.0
      vhem = vhem - vin
      vhem = vhem**2
      prob = 0.0
      probbreak = 0.0
      MS = 0
      C1 = m
      C3 = bint + 1
      C4 = bint + 2


C Calculate probability array and master connection matrix
C For master matrix MS, given B allowed binders, first B columns are possible 
C binding pairs. Second to last column is number of initial binders, last column
C is used to keep track of attempts and is initially equal to current binders
C Prob is probability of forming bond based on volume overlap
C probbreak is probability of breaking for each bound pair based on 
C force dependent off rate 
      do I=1, m
         array_1(I) = I

        do J=I+1,m
           prob(I,J) = 1 - exp(- kon* vollap(I,J)/vhem)
           probbreak(I,J) = 1 - exp(-fbreak(I,J))
           
           if (rmx(I,J) .GT. 0.5) then
                A1 = MS(I,C3)+1
                A2 = MS(J,C3)+1
                MS(I,A1) = J
                MS(J,A2) = I
                MS(I,C3) = A1
                MS(J,C3) = A2
                MS(I,C4) = A1
                MS(J,C4) = A2
           end if 
        end do 
      end do
 
        prob = (prob + TRANSPOSE(prob))
        probbreak = (probbreak+TRANSPOSE(probbreak))
C Randomly samples through each array until all beads have been sampled
C or max number of binders are formed 
C Bond formation 
      do while (C1 .GT. 0)

        CALL Random_Number(r3)
        r1 = floor(C1*r3) +1
        B1 = array_1(r1)
    
        array_2 = array_1
        C2 = C1 
            
           do while ((MS(B1,C4) .LT. bint) .AND. (C2 .GT. 0) )
 
                CALL Random_Number(r3)
                r2 = floor (C2*r3)+1
                B2 = array_2(r2)
C Condition checks that beads are greater than 1 bead away, the second bead 
C has not maxed its binding pairs, and they are not already bound 
              if ((ABS(B1 - B2) .GT. 1) .AND. (MS(B2,C4) .LT. bint)
     +  .AND. (rmxout(B1,B2) .LT. 0.5 )  ) then

                     CALL Random_Number(r3)
                        if( r3 .LT. prob(B1,B2)) then 
                                
                           MS(B1,C4) = MS(B1,C4)+1
                           MS(B2,C4) = MS(B2,C4)+1

                            rmxout(B1, B2) = 1
                            rmxout(B2, B1) = 1                          
                        end if 
              end if 
                 
              array_2(r2) = array_2(C2)
              C2 = C2 - 1

           end do 
C Bond breaking           
            MS(B1,C4) = bint
          A4 = MS(B1,C3)
          do I=1,A4
            B2 = MS(B1,I)
            if (MS(B2,1) . GT. 0) then
               CALL Random_Number(r3)
               if (r3 .LT. probbreak(B1,B2)) then 
                    rmxout(B1, B2) = 0
                    rmxout(B2, B1) = 0
               end if 
            end if 
           MS(B1,I) = 0
          end do
C move elements in array so B1 is not sampled again        
       
        array_1(r1) = array_1(C1)
        C1 = C1 - 1
      end do
      RETURN    
      END