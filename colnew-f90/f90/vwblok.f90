!*==VWBLOK.spg  processed by SPAG 6.72Dc at 06:36 on  8 Sep 2021
      SUBROUTINE VWBLOK(Xcol,Hrho,Jj,Wi,Vi,Ipvtw,Kd,Zval,Df,Acol,Dmzo,  &
                      & Ncomp,DFSUB,Msing)
!
!**********************************************************************
!
!c VWBLOK constructs a group of NCOMP rows of the WI and VI matrices.
!
!   purpose:
!
!      construct a group of  ncomp  rows of the matrices  wi  and  vi.
!      corresponding to an interior collocation point.
!
!
!   variables:
!
!      xcol   - the location of the collocation point.
!      jj     - xcol is the jj-th of k collocation points
!               in the i-th subinterval.
!      wi,vi  - the i-th block of the collocation matrix
!               before parameter condensation.
!      kd     - no. of rows in vi and wi .
!      zval   - z(xcol)
!      df     - the jacobian at xcol .
!      jcomp  - counter for the component being dealt with.
!
!**********************************************************************
      IMPLICIT NONE
!*--VWBLOK29
!*** Start of declarations inserted by SPAG
      REAL*8 Acol , ajl , basm , bl , Df , Dmzo , fact , ha , Hrho ,    &
           & paramsh , paramsl , Vi , Wi , Xcol , Zval
      INTEGER i0 , i1 , i2 , ICAre , id , IGUess , Ipvtw , ir , ITEr ,  &
            & iw , j , jcol , jcomp , jdf , Jj , jn , jv , jw , K , Kd
      INTEGER KDUm , l , LIMit , ll , lp1 , M , mj , MMAx , Msing ,     &
            & MSTar , NCDum , Ncomp , NONlin
!*** End of declarations inserted by SPAG
      DIMENSION Wi(Kd,1) , Vi(Kd,1) , Zval(1) , Dmzo(1) , Df(Ncomp,1)
      DIMENSION Ipvtw(1) , ha(7,4) , Acol(7,4) , basm(5)
!
      COMMON /COLORD/ K , NCDum , MSTar , KDUm , MMAx , M(20)
      COMMON /COLNLN/ NONlin , ITEr , LIMit , ICAre , IGUess
!
!...  if jj = 1 initialize  wi .
!
      IF ( Jj.LE.1 ) THEN
         DO id = 1 , Kd
            Wi(id,id) = 1.D0
         ENDDO
      ENDIF
!
!...  calculate local basis
!
      fact = 1.D0
      DO l = 1 , MMAx
         fact = fact*Hrho/DFLOAT(l)
         basm(l) = fact
         DO j = 1 , K
            ha(j,l) = fact*Acol(j,l)
         ENDDO
      ENDDO
!
!... zero jacobian
!
      DO jcol = 1 , MSTar
         DO ir = 1 , Ncomp
            Df(ir,jcol) = 0.D0
         ENDDO
      ENDDO
!
!...  build ncomp rows for interior collocation point x.
!...  the linear expressions to be constructed are:
!...   (m(id))
!...  u     -  df(id,1)*z(1) - ... - df(id,mstar)*z(mstar)
!...   id
!...  for id = 1 to ncomp.
!
      CALL DFSUB(paramsl,paramsh,Xcol,Zval,Df)
!      type (params) :: paramsL,paramsH
 
      i0 = (Jj-1)*Ncomp
      i1 = i0 + 1
      i2 = i0 + Ncomp
!
!...  evaluate  dmzo = dmz - df * zval  once for a new mesh
!
      IF ( NONlin.NE.0 .AND. ITEr.LE.0 ) THEN
         DO j = 1 , MSTar
            fact = -Zval(j)
            DO id = 1 , Ncomp
               Dmzo(i0+id) = Dmzo(i0+id) + fact*Df(id,j)
            ENDDO
         ENDDO
      ENDIF
!
!...  loop over the  ncomp  expressions to be set up for the
!...  current collocation point.
!
      DO j = 1 , MSTar
         DO id = 1 , Ncomp
            Vi(i0+id,j) = Df(id,j)
         ENDDO
      ENDDO
      jn = 1
      DO jcomp = 1 , Ncomp
         mj = M(jcomp)
         jn = jn + mj
         DO l = 1 , mj
            jv = jn - l
            jw = jcomp
            DO j = 1 , K
               ajl = -ha(j,l)
               DO iw = i1 , i2
                  Wi(iw,jw) = Wi(iw,jw) + ajl*Vi(iw,jv)
               ENDDO
               jw = jw + Ncomp
            ENDDO
            lp1 = l + 1
            IF ( l.NE.mj ) THEN
               DO ll = lp1 , mj
                  jdf = jn - ll
                  bl = basm(ll-l)
                  DO iw = i1 , i2
                     Vi(iw,jv) = Vi(iw,jv) + bl*Vi(iw,jdf)
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      IF ( Jj.LT.K ) RETURN
!
!   ...decompose the wi block and solve for the mstar columns of vi
!
!
!...  do parameter condensation
!
      Msing = 0
      CALL DGEFA(Wi,Kd,Kd,Ipvtw,Msing)
!
!...   check for singularity
!
      IF ( Msing.NE.0 ) RETURN
      DO j = 1 , MSTar
         CALL DGESL(Wi,Kd,Kd,Ipvtw,Vi(1,j),0)
      ENDDO
      END
