!*==GBLOCK.spg  processed by SPAG 6.72Dc at 07:09 on  8 Sep 2021
      SUBROUTINE GBLOCK(H,Gi,Nrow,Irow,Wi,Vi,Kd,Rhsz,Rhsdmz,Ipvtw,Mode)
!
!**********************************************************************
!
!c GBLOCK constructs certain rows of the collocation matrix.
!
!   purpose:
!
!      construct collocation matrix rows according to mode:
!      mode = 1  -  a group of  mstar    rows corresponding
!                   an interior mesh interval.
!           = 2  -  corresponding right hand side
!
!   variables:
!
!      h      - the  local stepsize.
!      gi     - the sub-block of the collocation matrix in
!               which the equations are to be formed.
!      wi     - the sub-block of noncondensed collocation equations,
!               left-hand side part.
!      vi     - the sub-block of noncondensed collocation equations,
!               right-hand side part.
!      rhsdmz - the inhomogenous term of the uncondensed collocation
!               equations.
!      rhsz   - the inhomogenous term of the condensed collocation
!               equations.
!      nrow   - no. of rows in gi.
!      irow   - the first row in gi to be used for equations.
!
!**********************************************************************
      IMPLICIT NONE
!*--GBLOCK33
!*** Start of declarations inserted by SPAG
      REAL*8 ACOl , ASAve , B , basm , fact , Gi , H , hb , Rhsdmz ,    &
           & Rhsz , rsum , Vi , Wi
      INTEGER icomp , id , ind , Ipvtw , ir , Irow , j , jcol , jcomp , &
            & jd , K , Kd , KDUm , l , ll , M , mj , MMAx , Mode , MSTar
      INTEGER NCOmp , Nrow
!*** End of declarations inserted by SPAG
      DIMENSION hb(7,4) , basm(5)
      DIMENSION Gi(Nrow,1) , Wi(1) , Vi(Kd,1)
      DIMENSION Rhsz(1) , Rhsdmz(1) , Ipvtw(1)
!
      COMMON /COLORD/ K , NCOmp , MSTar , KDUm , MMAx , M(20)
      COMMON /COLBAS/ B(7,4) , ACOl(28,7) , ASAve(28,4)
!
!...  compute local basis
!
      fact = 1.D0
      basm(1) = 1.D0
      DO l = 1 , MMAx
         fact = fact*H/DFLOAT(l)
         basm(l+1) = fact
         DO j = 1 , K
            hb(j,l) = fact*B(j,l)
         ENDDO
      ENDDO
!
!...  branch according to  m o d e
!
      IF ( Mode.EQ.2 ) THEN
!
!...  compute the appropriate piece of  rhsz
!
         CALL DGESL(Wi,Kd,Kd,Ipvtw,Rhsdmz,0)
         ir = Irow
         DO jcomp = 1 , NCOmp
            mj = M(jcomp)
            ir = ir + mj
            DO l = 1 , mj
               ind = jcomp
               rsum = 0.D0
               DO j = 1 , K
                  rsum = rsum + hb(j,l)*Rhsdmz(ind)
                  ind = ind + NCOmp
               ENDDO
               Rhsz(ir-l) = rsum
            ENDDO
         ENDDO
         GOTO 99999
      ENDIF
!
!...  set right gi-block to identity
!
      DO j = 1 , MSTar
         DO ir = 1 , MSTar
            Gi(Irow-1+ir,j) = 0.D0
            Gi(Irow-1+ir,MSTar+j) = 0.D0
         ENDDO
         Gi(Irow-1+j,MSTar+j) = 1.D0
      ENDDO
!
!...  compute the block gi
!
      ir = Irow
      DO icomp = 1 , NCOmp
         mj = M(icomp)
         ir = ir + mj
         DO l = 1 , mj
            id = ir - l
            DO jcol = 1 , MSTar
               ind = icomp
               rsum = 0.D0
               DO j = 1 , K
                  rsum = rsum - hb(j,l)*Vi(ind,jcol)
                  ind = ind + NCOmp
               ENDDO
               Gi(id,jcol) = rsum
            ENDDO
            jd = id - Irow
            DO ll = 1 , l
               Gi(id,jd+ll) = Gi(id,jd+ll) - basm(ll)
            ENDDO
         ENDDO
      ENDDO
      RETURN
99999 END
