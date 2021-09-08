!*==ERRCHK.spg  processed by SPAG 6.72Dc at 07:06 on  8 Sep 2021
      SUBROUTINE ERRCHK(Xi,Z,Dmz,Valstr,Ifin)
!
!**********************************************************************
!
!c ERRCHK determines error estimates and tests error tolerances.
!
!      purpose
!               determine the error estimates and test to see if the
!               error tolerances are satisfied.
!
!      variables
!        xi     - current mesh points
!        valstr - values of the previous solution which are needed
!                 for the extrapolation- like error estimate.
!        wgterr - weights used in the extrapolation-like error
!                 estimate. the array values are assigned in
!                 subroutine  consts.
!        errest - storage for error estimates
!        err    - temporary storage used for error estimates
!        z      - approximate solution on mesh xi
!        ifin   - a 0-1 variable. on return it indicates whether
!                 the error tolerances were satisfied
!        mshflg - is set by errchk to indicate to newmsh whether
!                 any values of the current solution are stored in
!                 the array valstr. (0 for no, 1 for yes)
!
!**********************************************************************
!
      IMPLICIT NONE
!*--ERRCHK31
!*** Start of declarations inserted by SPAG
      REAL*8 ACOl , ASAve , B , Dmz , dummy , err , errest , PREcis ,   &
           & ROOt , TOL , TOLin , Valstr , WGTerr , WGTmsh , x , Xi , Z
      INTEGER i , iback , Ifin , IOUt , IPRint , j , JTOl , K , KD ,    &
            & knew , kstore , l , lj , ltjz , LTOl , ltolj , M , mj ,   &
            & MMAx , MSHalt
      INTEGER MSHflg , MSHlmt , MSHnum , MSTar , N , NCOmp , NDMz ,     &
            & NMAx , NOLd , NTOl , NZ
!*** End of declarations inserted by SPAG
      DIMENSION err(40) , errest(40) , dummy(1)
      DIMENSION Xi(1) , Z(1) , Dmz(1) , Valstr(1)
!
      COMMON /COLOUT/ PREcis , IOUt , IPRint
      COMMON /COLORD/ K , NCOmp , MSTar , KD , MMAx , M(20)
      COMMON /COLAPR/ N , NOLd , NMAx , NZ , NDMz
      COMMON /COLMSH/ MSHflg , MSHnum , MSHlmt , MSHalt
      COMMON /COLBAS/ B(28) , ACOl(28,7) , ASAve(28,4)
      COMMON /COLEST/ TOL(40) , WGTmsh(40) , WGTerr(40) , TOLin(40) ,   &
                    & ROOt(40) , JTOl(40) , LTOl(40) , NTOl
!
!...  error estimates are to be generated and tested
!...  to see if the tolerance requirements are satisfied.
!
      Ifin = 1
      MSHflg = 1
      DO j = 1 , MSTar
         errest(j) = 0.D0
      ENDDO
      DO iback = 1 , N
         i = N + 1 - iback
!
!...       the error estimates are obtained by combining values of
!...       the numerical solutions for two meshes.
!...       for each value of iback we will consider the two
!...       approximations at 2 points in each of
!...       the new subintervals.  we work backwards through
!...       the subinterval so that new values can be stored
!...       in valstr in case they prove to be needed later
!...       for an error estimate. the routine  newmsh
!...       filled in the needed values of the old solution
!...       in valstr.
!
         knew = (4*(i-1)+2)*MSTar + 1
         kstore = (2*(i-1)+1)*MSTar + 1
         x = Xi(i) + (Xi(i+1)-Xi(i))*2.D0/3.D0
         CALL APPROX(i,x,Valstr(knew),ASAve(1,3),dummy,Xi,N,Z,Dmz,K,    &
                   & NCOmp,MMAx,M,MSTar,4,dummy,0)
         DO l = 1 , MSTar
            err(l) = WGTerr(l)*DABS(Valstr(knew)-Valstr(kstore))
            knew = knew + 1
            kstore = kstore + 1
         ENDDO
         knew = (4*(i-1)+1)*MSTar + 1
         kstore = 2*(i-1)*MSTar + 1
         x = Xi(i) + (Xi(i+1)-Xi(i))/3.D0
         CALL APPROX(i,x,Valstr(knew),ASAve(1,2),dummy,Xi,N,Z,Dmz,K,    &
                   & NCOmp,MMAx,M,MSTar,4,dummy,0)
         DO l = 1 , MSTar
            err(l) = err(l) + WGTerr(l)                                 &
                   & *DABS(Valstr(knew)-Valstr(kstore))
            knew = knew + 1
            kstore = kstore + 1
         ENDDO
!
!...       find component-wise maximum error estimate
!
         DO l = 1 , MSTar
            errest(l) = DMAX1(errest(l),err(l))
         ENDDO
!
 
 
!...       test whether the tolerance requirements are satisfied
!...       in the i-th interval.
!
         IF ( Ifin.NE.0 ) THEN
            DO j = 1 , NTOl
               ltolj = LTOl(j)
               ltjz = ltolj + (i-1)*MSTar
               IF ( err(ltolj).GT.TOLin(j)*(DABS(Z(ltjz))+1.D0) )       &
                  & Ifin = 0
            ENDDO
         ENDIF
      ENDDO
      IF ( IPRint.GE.0 ) RETURN
      WRITE (IOUt,99001)
99001 FORMAT (/' THE ESTIMATED ERRORS ARE,')
      lj = 1
      DO j = 1 , NCOmp
         mj = lj - 1 + M(j)
         WRITE (IOUt,99002) j , (errest(l),l=lj,mj)
!--------------------------------------------------------------
99002    FORMAT (' U(',I2,') -',4D12.4)
         lj = mj + 1
      ENDDO
      RETURN
      END
