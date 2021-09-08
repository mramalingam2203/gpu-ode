!*==CONSTS.spg  processed by SPAG 6.72Dc at 05:10 on  8 Sep 2021
      SUBROUTINE CONSTS(K,Rho,Coef)
!
!**********************************************************************
!
!c CONSTS assigns values to various array constants.
!
!   purpose
!            assign (once) values to various array constants.
!
!   arrays assigned during compilation:
!     cnsts1 - weights for extrapolation error estimate
!     cnsts2 - weights for mesh selection
!              (the above weights come from the theoretical form for
!              the collocation error -- see [3])
!
!   arrays assigned during execution:
!     wgterr - the particular values of cnsts1 used for current run
!              (depending on k, m)
!     wgtmsh - gotten from the values of cnsts2 which in turn are
!              the constants in the theoretical expression for the
!              errors. the quantities in wgtmsh are 10x the values
!              in cnsts2 so that the mesh selection algorithm
!              is aiming for errors .1x as large as the user
!              requested tolerances.
!     jtol   - components of differential system to which tolerances
!              refer (viz, if ltol(i) refers to a derivative of u(j),
!              then jtol(i)=j)
!     root   - reciprocals of expected rates of convergence of compo-
!              nents of z(j) for which tolerances are specified
!     rho    - the k collocation points on (0,1)
!     coef   -
!     acol  -  the runge-kutta coefficients values at collocation
!              points
!
!**********************************************************************
!
      IMPLICIT NONE
!*--CONSTS39
!*** Start of declarations inserted by SPAG
      REAL*8 ACOl , ASAve , B , cnsts1 , cnsts2 , Coef , dummy , Rho ,  &
     &       ROOt , TOL , TOLin , WGTerr , WGTmsh
      INTEGER i , iz , j , jcomp , JTOl , K , KD , KDUm , koff , l ,    &
     &        LTOl , ltoli , M , mj , MMAx , MSTar , mtot , NCOmp , NTOl
!*** End of declarations inserted by SPAG
      DIMENSION Rho(7) , Coef(K,1) , cnsts1(28) , cnsts2(28) , dummy(1)
!
      COMMON /COLORD/ KDUm , NCOmp , MSTar , KD , MMAx , M(20)
      COMMON /COLBAS/ B(28) , ACOl(28,7) , ASAve(28,4)
      COMMON /COLEST/ TOL(40) , WGTmsh(40) , WGTerr(40) , TOLin(40) ,   &
     &                ROOt(40) , JTOl(40) , LTOl(40) , NTOl
!
      DATA cnsts1/.25D0 , .625D-1 , 7.2169D-2 , 1.8342D-2 , 1.9065D-2 , &
     &     5.8190D-2 , 5.4658D-3 , 5.3370D-3 , 1.8890D-2 , 2.7792D-2 ,  &
     &     1.6095D-3 , 1.4964D-3 , 7.5938D-3 , 5.7573D-3 , 1.8342D-2 ,  &
     &     4.673D-3 , 4.150D-4 , 1.919D-3 , 1.468D-3 , 6.371D-3 ,       &
     &     4.610D-3 , 1.342D-4 , 1.138D-4 , 4.889D-4 , 4.177D-4 ,       &
     &     1.374D-3 , 1.654D-3 , 2.863D-3/
      DATA cnsts2/1.25D-1 , 2.604D-3 , 8.019D-3 , 2.170D-5 , 7.453D-5 , &
     &     5.208D-4 , 9.689D-8 , 3.689D-7 , 3.100D-6 , 2.451D-5 ,       &
     &     2.691D-10 , 1.120D-9 , 1.076D-8 , 9.405D-8 , 1.033D-6 ,      &
     &     5.097D-13 , 2.290D-12 , 2.446D-11 , 2.331D-10 , 2.936D-9 ,   &
     &     3.593D-8 , 7.001D-16 , 3.363D-15 , 3.921D-14 , 4.028D-13 ,   &
     &     5.646D-12 , 7.531D-11 , 1.129D-9/
!
!...  assign weights for error estimate
!
      koff = K*(K+1)/2
      iz = 1
      DO j = 1 , NCOmp
         mj = M(j)
         DO l = 1 , mj
            WGTerr(iz) = cnsts1(koff-mj+l)
            iz = iz + 1
         ENDDO
      ENDDO
!
!...  assign array values for mesh selection: wgtmsh, jtol, and root
!
      jcomp = 1
      mtot = M(1)
      DO i = 1 , NTOl
         ltoli = LTOl(i)
 50      IF ( ltoli.LE.mtot ) THEN
            JTOl(i) = jcomp
            WGTmsh(i) = 1.D1*cnsts2(koff+ltoli-mtot)/TOLin(i)
            ROOt(i) = 1.D0/DFLOAT(K+mtot-ltoli+1)
         ELSE
            jcomp = jcomp + 1
            mtot = mtot + M(jcomp)
            GOTO 50
         ENDIF
      ENDDO
!
!...  specify collocation points
!
      IF ( K.EQ.2 ) THEN
         Rho(2) = .57735026918962576451D0
         Rho(1) = -Rho(2)
      ELSEIF ( K.EQ.3 ) THEN
         Rho(3) = .77459666924148337704D0
         Rho(2) = .0D0
         Rho(1) = -Rho(3)
      ELSEIF ( K.EQ.4 ) THEN
         Rho(4) = .86113631159405257523D0
         Rho(3) = .33998104358485626480D0
         Rho(2) = -Rho(3)
         Rho(1) = -Rho(4)
      ELSEIF ( K.EQ.5 ) THEN
         Rho(5) = .90617984593866399280D0
         Rho(4) = .53846931010568309104D0
         Rho(3) = .0D0
         Rho(2) = -Rho(4)
         Rho(1) = -Rho(5)
      ELSEIF ( K.EQ.6 ) THEN
         Rho(6) = .93246951420315202781D0
         Rho(5) = .66120938646626451366D0
         Rho(4) = .23861918608319690863D0
         Rho(3) = -Rho(4)
         Rho(2) = -Rho(5)
         Rho(1) = -Rho(6)
      ELSEIF ( K.EQ.7 ) THEN
         Rho(7) = .949107991234275852452D0
         Rho(6) = .74153118559939443986D0
         Rho(5) = .40584515137739716690D0
         Rho(4) = 0.D0
         Rho(3) = -Rho(5)
         Rho(2) = -Rho(6)
         Rho(1) = -Rho(7)
      ELSE
         Rho(1) = 0.D0
      ENDIF
!
!...  map (-1,1) to (0,1) by  t = .5 * (1. + x)
!
      DO j = 1 , K
         Rho(j) = .5D0*(1.D0+Rho(j))
      ENDDO
!
!...  now find runge-kutta coeffitients b, acol and asave
!...  the values of asave are to be used in  newmsh  and errchk .
!
      DO j = 1 , K
         DO i = 1 , K
            Coef(i,j) = 0.D0
         ENDDO
         Coef(j,j) = 1.D0
         CALL VMONDE(Rho,Coef(1,j),K)
      ENDDO
      CALL RKBAS(1.D0,Coef,K,MMAx,B,dummy,0)
      DO i = 1 , K
         CALL RKBAS(Rho(i),Coef,K,MMAx,ACOl(1,i),dummy,0)
      ENDDO
      CALL RKBAS(1.D0/6.D0,Coef,K,MMAx,ASAve(1,1),dummy,0)
      CALL RKBAS(1.D0/3.D0,Coef,K,MMAx,ASAve(1,2),dummy,0)
      CALL RKBAS(2.D0/3.D0,Coef,K,MMAx,ASAve(1,3),dummy,0)
      CALL RKBAS(5.D0/6.D0,Coef,K,MMAx,ASAve(1,4),dummy,0)
      END
