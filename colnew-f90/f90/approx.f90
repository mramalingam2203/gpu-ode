!*==APPROX.spg  processed by SPAG 6.72Dc at 05:08 on  8 Sep 2021
      SUBROUTINE APPROX(I,X,Zval,A,Coef,Xi,N,Z,Dmz,K,Ncomp,Mmax,M,Mstar,&
     &                  Mode,Dmval,Modm)
!
!**********************************************************************
!
!c APPROX evaluates the computes solution at a point.
!
!   purpose
!                                    (1)       (m1-1)     (mncomp-1)
!           evaluate z(u(x))=(u (x),u (x),...,u  (x),...,u  (x)      )
!                              1     1         1          mncomp
!           at one point x.
!
!   variables
!     a      - array of mesh independent rk-basis coefficients
!     basm   - array of mesh dependent monomial coefficients
!     xi     - the current mesh (having n subintervals)
!     z      - the current solution vector
!     dmz    - the array of mj-th derivatives of the current solution
!     mode   - determines the amount of initialization needed
!            = 4  forms z(u(x)) using z, dmz and ha
!            = 3  as in =4, but computes local rk-basis
!            = 2  as in =3, but determines i such that
!                       xi(i) .le. x .lt. xi(i+1) (unless x=xi(n+1))
!            = 1  retrieve  z=z(u(x(i)))  directly
!
!**********************************************************************
!
      IMPLICIT NONE
!*--APPROX31
!*** Start of declarations inserted by SPAG
      REAL*8 A , bm , Coef , dm , Dmval , Dmz , fact , PREcis , s , X , &
     &       Xi , Z , zsum , Zval
      INTEGER I , idmz , ileft , ind , IOUt , IPRint , ir , iright ,    &
     &        iz , j , jcomp , K , l , lb , ll , M , mj , Mmax , Mode , &
     &        Modm
      INTEGER Mstar , N , Ncomp
!*** End of declarations inserted by SPAG
      DIMENSION Zval(1) , Dmval(1) , Xi(1) , M(1) , A(7,1) , dm(7)
      DIMENSION Z(1) , Dmz(1) , bm(4) , Coef(1)
!
      COMMON /COLOUT/ PREcis , IOUt , IPRint
!
      IF ( Mode.EQ.2 ) THEN
!
!...  mode = 2 ,  locate i so  xi(i) .le. x .lt. xi(i+1)
!
         IF ( X.LT.Xi(1)-PREcis .OR. X.GT.Xi(N+1)+PREcis ) THEN
            IF ( IPRint.LT.1 ) WRITE (IOUt,99001) X , Xi(1) , Xi(N+1)
!--------------------------------------------------------------------
99001       FORMAT (' ****** DOMAIN ERROR IN APPROX ******'/' X =',     &
     &              D20.10,'   ALEFT =',D20.10,'   ARIGHT =',D20.10)
            IF ( X.LT.Xi(1) ) X = Xi(1)
            IF ( X.GT.Xi(N+1) ) X = Xi(N+1)
         ENDIF
         IF ( I.GT.N .OR. I.LT.1 ) I = (N+1)/2
         ileft = I
         IF ( X.LT.Xi(ileft) ) THEN
            iright = ileft - 1
            DO l = 1 , iright
               I = iright + 1 - l
               IF ( X.GE.Xi(I) ) GOTO 100
            ENDDO
         ELSE
            DO l = ileft , N
               I = l
               IF ( X.LT.Xi(l+1) ) GOTO 100
            ENDDO
         ENDIF
      ELSEIF ( Mode.EQ.3 ) THEN
      ELSEIF ( Mode.EQ.4 ) THEN
         GOTO 200
      ELSE
!
!...  mode = 1 , retrieve  z( u(x) )  directly for x = xi(i).
!
         X = Xi(I)
         iz = (I-1)*Mstar
         DO j = 1 , Mstar
            iz = iz + 1
            Zval(j) = Z(iz)
         ENDDO
         RETURN
      ENDIF
!
!...  mode = 2 or 3 , compute mesh independent rk-basis.
!
 100  s = (X-Xi(I))/(Xi(I+1)-Xi(I))
      CALL RKBAS(s,Coef,K,Mmax,A,dm,Modm)
!
!...  mode = 2, 3, or 4 , compute mesh dependent rk-basis.
!
 200  bm(1) = X - Xi(I)
      DO l = 2 , Mmax
         bm(l) = bm(1)/DFLOAT(l)
      ENDDO
!
!...  evaluate  z( u(x) ).
!
      ir = 1
      iz = (I-1)*Mstar + 1
      idmz = (I-1)*K*Ncomp
      DO jcomp = 1 , Ncomp
         mj = M(jcomp)
         ir = ir + mj
         iz = iz + mj
         DO l = 1 , mj
            ind = idmz + jcomp
            zsum = 0.D0
            DO j = 1 , K
               zsum = zsum + A(j,l)*Dmz(ind)
               ind = ind + Ncomp
            ENDDO
            DO ll = 1 , l
               lb = l + 1 - ll
               zsum = zsum*bm(lb) + Z(iz-ll)
            ENDDO
            Zval(ir-l) = zsum
         ENDDO
      ENDDO
      IF ( Modm.EQ.0 ) RETURN
!
!...  for modm = 1 evaluate  dmval(j) = mj-th derivative of uj.
!
      DO jcomp = 1 , Ncomp
         Dmval(jcomp) = 0.D0
      ENDDO
      idmz = idmz + 1
      DO j = 1 , K
         fact = dm(j)
         DO jcomp = 1 , Ncomp
            Dmval(jcomp) = Dmval(jcomp) + fact*Dmz(idmz)
            idmz = idmz + 1
         ENDDO
      ENDDO
      RETURN
      END
