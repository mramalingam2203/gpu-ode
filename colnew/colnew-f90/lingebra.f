
subroutine vmonde ( rho, coef, k )

c
c**********************************************************************
c
cc vmonde solves a vandermonde linear system.
c
c   purpose
c          solve vandermonde system v * x = e
c          with  v(i,j) = rho(j)**(i-1)/(i-1)! .
c
c**********************************************************************
c
      integer k, i,ifac,j,km1,kmi
      double precision rho(k), coef(k)
c
      if ( k .eq. 1 )                             return
      km1 = k - 1
      do 10 i = 1, km1
         kmi = k - i
         do 10 j = 1, kmi
           coef(j) = (coef(j+1) - coef(j)) / (rho(j+i) - rho(j))
  10  continue
c
      ifac = 1
      do 40 i = 1, km1
         kmi = k + 1 - i
         do 30 j = 2, kmi
  30        coef(j) = coef(j) - rho(j+i-1) * coef(j-1)
         coef(kmi) = dfloat(ifac) * coef(kmi)
         ifac = ifac * i
  40  continue
      coef(1) = dfloat(ifac) * coef(1)
      return
      end
      subroutine horder (i, uhigh, hi, dmz, ncomp, k)
c
c**********************************************************************
c
cc horder determines the highest order derivatives of the current solution.
c
c   purpose
c           determine highest order (piecewise constant) derivatives
c           of the current collocation solution
c
c   variables
c     hi     - the stepsize, hi = xi(i+1) - xi(i)
c     dmz    - vector of mj-th derivative of the solution
c     uhigh  - the array of highest order (piecewise constant)
c              derivatives of the approximate solution on
c              (xi(i),xi(i+1)), viz,
c                          (k+mj-1)
c              uhigh(j) = u   (x)    on (xi(i),xi(i+1))
c                          j
c
c**********************************************************************
c
      implicit real*8 (a-h,o-z)
      dimension uhigh(1), dmz(1)
c
      common /colloc/ rho(7), coef(49)
c
      dn = 1.d0 / hi**(k-1)
c
c...  loop over the ncomp solution components
c
      do 10 id = 1, ncomp
         uhigh(id) = 0.d0
   10 continue
      kin = 1
      idmz = (i-1) * k * ncomp + 1
      do 30 j = 1, k
         fact = dn * coef(kin)
         do 20 id = 1, ncomp
            uhigh(id) = uhigh(id)  +  fact * dmz(idmz)
            idmz = idmz + 1
   20    continue
         kin = kin + k
   30 continue
      return
      end
      subroutine dmzsol (kd, mstar, n, v, z, dmz)
c
c**********************************************************************
c
cc dmzsol computes dmz in a blockwise manner.
c
c   purpose
c          compute dmz in a blockwise manner
c          dmz(i) = dmz(i)  +  v(i) * z(i), i = 1,...,n
c
c**********************************************************************
c
      implicit real*8 (a-h,o-z)
      dimension v(kd,1), dmz(kd,1), z(1)
c
      jz = 1
      do 30 i = 1, n
         do 20 j = 1, mstar
            fact = z(jz)
            do 10 l = 1, kd
               dmz(l,i) = dmz(l,i)  +  fact * v(l,jz)
   10       continue
            jz = jz + 1
   20    continue
   30 continue
      return
      end
c----------------------------------------------------------------------
c                            p a r t  5
c          we list here a modified (column oriented, faster)
c          version of the package solveblok of de boor - weiss [5].
c          we also give a listing of the linpack
c          routines dgefa und dgesl used by colnew.
c----------------------------------------------------------------------
c
      subroutine fcblok (bloks, integs, nbloks, ipivot, scrtch, info)

c**********************************************************************
c
cc fcblok calls subroutines factrb and shiftb.
c
c     fcblok  supervises the plu factorization with pivoting of
c     scaled rows of the almost block diagonal matrix stored in the
c     arrays  bloks  and  integs .
c
c     factrb = subprogram which carries out steps 1,...,last of gauss
c            elimination (with pivoting) for an individual block.
c     shiftb = subprogram which shifts the remaining rows to the top of
c            the next block
c
c     parameters
c      bloks   an array that initially contains the almost block diago-
c            nal matrix  a  to be factored, and on return contains the
c            computed factorization of  a .
c      integs  an integer array describing the block structure of  a .
c      nbloks  the number of blocks in  a .
c      ipivot  an integer array of dimension   sum (integs(3,n) ; n=1,
c            ...,nbloks) which, on return, contains the pivoting stra-
c            tegy used.
c      scrtch  work area required, of length  max (integs(1,n) ; n=1,
c            ...,nbloks).
c      info    output parameter;
c            = 0  in case matrix was found to be nonsingular.
c            otherwise,
c            = n if the pivot element in the nth gauss step is zero.
c
c**********************************************************************
c
      integer integs(3,nbloks),ipivot(1),info, i,index,indexn,last,
     1        ncol,nrow
      double precision bloks(1),scrtch(1)
      info = 0
      indexx = 1
      indexn = 1
      i = 1
c
c...  loop over the blocks.  i  is loop index
c
   10      index = indexn
           nrow = integs(1,i)
           ncol = integs(2,i)
           last = integs(3,i)
c
c...       carry out elimination on the i-th block until next block
c...       enters, i.e., for columns 1,...,last  of i-th block.
c
           call factrb ( bloks(index), ipivot(indexx), scrtch, nrow,
     1                   ncol, last, info)
c
c...       check for having reached a singular block or the last block
c
           if ( info .ne. 0 )                       go to 20
           if ( i .eq. nbloks )                     return
           i = i+1
           indexn = nrow * ncol + index
           indexx = indexx + last
c
c...       put the rest of the i-th block onto the next block
c
           call shiftb ( bloks(index), nrow, ncol, last,
     1                   bloks(indexn), integs(1,i), integs(2,i) )
      go to 10
   20 info = info + indexx - 1
      return
      end
      subroutine factrb ( w, ipivot, d, nrow, ncol, last, info)
c
c********************************************************************
c
cc factrb constructs a partial plu factorization of a matrix.
c
c     adapted from p.132 of  element.numer.analysis  by conte-de boor
c
c     constructs a partial plu factorization, corresponding to steps
c      1,..., last   in gauss elimination, for the matrix  w  of
c      order ( nrow ,  ncol ), using pivoting of scaled rows.
c
c     parameters
c       w       contains the (nrow,ncol) matrix to be partially factored
c               on input, and the partial factorization on output.
c       ipivot  an integer array of length last containing a record of
c               the pivoting strategy used; explicit interchanges
c               are used for pivoting.
c       d       a work array of length nrow used to store row sizes
c               temporarily.
c       nrow    number of rows of w.
c       ncol    number of columns of w.
c       last    number of elimination steps to be carried out.
c       info    on output, zero if the matrix is found to be non-
c               singular, in case a zero pivot was encountered in row
c               n,  info = n on output.
c
c**********************************************************************
c
      integer ipivot(nrow),ncol,last,info, i,j,k,l,kp1
      double precision w(nrow,ncol),d(nrow), colmax,t,s
      double precision dabs,dmax1
c
c...  initialize  d
c
      do 10 i = 1, nrow
        d(i) = 0.d0
   10 continue
      do 20 j = 1, ncol
        do 20 i = 1, nrow
          d(i) = dmax1( d(i) , dabs(w(i,j)))
   20 continue
c
c...  gauss elimination with pivoting of scaled rows, loop over
c...  k=1,.,last
c
      k = 1
c
c...  as pivot row for k-th step, pick among the rows not yet used,
c...  i.e., from rows  k ,..., nrow , the one whose k-th entry
c...  (compared to the row size) is largest. then, if this row
c...  does not turn out to be row k, interchange row k with this
c...  particular row and redefine ipivot(k).
c
   30      continue
           if ( d(k) .eq. 0.d0 )                    go to 90
           if (k .eq. nrow)                         go to 80
           l = k
           kp1 = k+1
           colmax = dabs(w(k,k)) / d(k)
c
c...       find the (relatively) largest pivot
c
           do 40 i = kp1, nrow
             if ( dabs(w(i,k)) .le. colmax * d(i) ) go to 40
             colmax = dabs(w(i,k)) / d(i)
             l = i
   40      continue
           ipivot(k) = l
           t = w(l,k)
           s = d(l)
           if ( l .eq. k )                          go to 50
             w(l,k) = w(k,k)
             w(k,k) = t
             d(l) = d(k)
             d(k) = s
   50      continue
c
c...       if pivot element is too small in absolute value, declare
c...       matrix to be noninvertible and quit.
c
           if ( dabs(t)+d(k) .le. d(k) )            go to 90
c
c...       otherwise, subtract the appropriate multiple of the pivot
c...       row from remaining rows, i.e., the rows (k+1),..., (nrow)
c...       to make k-th entry zero. save the multiplier in its place.
c...       for high performance do this operations column oriented.
c
           t = -1.0d0/t
           do 60 i = kp1, nrow
   60        w(i,k) = w(i,k) * t
           do 70 j=kp1,ncol
             t = w(l,j)
             if ( l .eq. k )                        go to 62
               w(l,j) = w(k,j)
               w(k,j) = t
   62        if ( t .eq. 0.d0 )                     go to 70
             do 64 i = kp1, nrow
   64           w(i,j) = w(i,j) + w(i,k) * t
   70      continue
           k = kp1
c
c...       check for having reached the next block.
c
           if ( k .le. last )                       go to 30
      return
c
c...  if  last  .eq. nrow , check now that pivot element in last row
c...  is nonzero.
c
   80 if( dabs(w(nrow,nrow))+d(nrow) .gt. d(nrow) ) return
c
c...  singularity flag set
c
   90 info = k
      return
      end
      subroutine shiftb (ai, nrowi, ncoli, last, ai1, nrowi1, ncoli1)
c
c*********************************************************************
c
cc shiftb shifts rows in the current block.
c
c     shifts the rows in current block, ai, not used as pivot rows, if
c     any, i.e., rows  (last+1),..., (nrowi), onto the first mmax =
c      = nrow-last  rows of the next block, ai1, with column last+j of
c      ai  going to column j , j=1,...,jmax=ncoli-last. the remaining
c     columns of these rows of ai1 are zeroed out.
c
c                                picture
c
c          original situation after         results in a new block i+1
c          last = 2 columns have been       created and ready to be
c          done in factrb (assuming no      factored by next factrb
c          interchanges of rows)            call.
c                      1
c                 x  x 1x  x  x           x  x  x  x  x
c                      1
c                 0  x 1x  x  x           0  x  x  x  x
c     block i          1                       ---------------
c     nrowi = 4   0  0 1x  x  x           0  0 1x  x  x  0  01
c     ncoli = 5        1                       1             1
c     last = 2    0  0 1x  x  x           0  0 1x  x  x  0  01
c     -------------------------------          1             1   new
c                      1x  x  x  x  x          1x  x  x  x  x1  block
c                      1                       1             1   i+1
c     block i+1        1x  x  x  x  x          1x  x  x  x  x1
c     nrowi1= 5        1                       1             1
c     ncoli1= 5        1x  x  x  x  x          1x  x  x  x  x1
c     -------------------------------          1-------------1
c                      1
c
c*********************************************************************
c
      integer last, j,jmax,jmaxp1,m,mmax
      double precision ai(nrowi,ncoli),ai1(nrowi1,ncoli1)
      mmax = nrowi - last
      jmax = ncoli - last
      if (mmax .lt. 1 .or. jmax .lt. 1)             return
c
c...  put the remainder of block i into ai1
c
      do 10 j=1,jmax
           do 10 m=1,mmax
   10 ai1(m,j) = ai(last+m,last+j)
      if (jmax .eq. ncoli1)                         return
c
c...  zero out the upper right corner of ai1
c
      jmaxp1 = jmax + 1
      do 20 j=jmaxp1,ncoli1
           do 20 m=1,mmax
   20 ai1(m,j) = 0.d0
      return
      end
      subroutine sbblok ( bloks, integs, nbloks, ipivot, x )
c
c**********************************************************************
c
cc sbblok calls subfor and subbak so solve a block linear system.
c
c     calls subroutines  subfor  and  subbak .
c
c     supervises the solution (by forward and backward substitution) of
c     the linear system  a*x = b  for x, with the plu factorization of
c     a  already generated in  fcblok .  individual blocks of
c     equations are solved via  subfor  and  subbak .
c
c    parameters
c       bloks, integs, nbloks, ipivot    are as on return from fcblok.
c       x       on input: the right hand side, in dense storage
c               on output: the solution vector
c
c*********************************************************************
c
      integer integs(3,nbloks),ipivot(1), i,index,indexx,j,last,
     1        nbp1,ncol,nrow
      double precision bloks(1), x(1)
c
c...  forward substitution pass
c
      index = 1
      indexx = 1
      do 10 i = 1, nbloks
           nrow = integs(1,i)
           last = integs(3,i)
           call subfor ( bloks(index), ipivot(indexx), nrow, last,
     1                   x(indexx) )
           index = nrow * integs(2,i) + index
   10 indexx = indexx + last
c
c...  back substitution pass
c
      nbp1 = nbloks + 1
      do 20 j = 1, nbloks
           i = nbp1 - j
           nrow = integs(1,i)
           ncol = integs(2,i)
           last = integs(3,i)
           index = index - nrow * ncol
           indexx = indexx - last
   20 call subbak ( bloks(index), nrow, ncol, last, x(indexx) )
      return
      end
      subroutine subfor ( w, ipivot, nrow, last, x )
c
c**********************************************************************
c
cc subfor carries out the forward pass of substitution for the current block.
c
c     carries out the forward pass of substitution for the current
c     block, i.e., the action on the right side corresponding to the
c     elimination carried out in  factrb  for this block.
c
c    parameters
c       w, ipivot, nrow, last  are as on return from factrb.
c       x(j)  is expected to contain, on input, the right side of j-th
c             equation for this block, j=1,...,nrow.
c       x(j)  contains, on output, the appropriately modified right
c             side of equation (j) in this block, j=1,...,last and
c             for j=last+1,...,nrow.
c
c*********************************************************************
c
      integer ipivot(last), ip,k,kp1,lstep
      double precision w(nrow,last), x(nrow), t
c
      if ( nrow .eq. 1 )                            return
      lstep = min0( nrow-1 , last )
      do 20 k = 1, lstep
           kp1 = k + 1
           ip = ipivot(k)
           t = x(ip)
           x(ip) = x(k)
           x(k) = t
           if ( t .eq. 0.d0 )                       go to 20
           do 10 i = kp1, nrow
   10         x(i) = x(i) + w(i,k) * t
   20 continue
   30 return
      end
      subroutine subbak ( w, nrow, ncol, last, x )
c
c*********************************************************************
c
cc subbak carries out backsubstitution for the current block.
c
c    parameters
c       w, ipivot, nrow, ncol, last  are as on return from factrb.
c       x(1),...,x(ncol)  contains, on input, the right side for the
c               equations in this block after backsubstitution has been
c               carried up to but not including equation (last).
c               means that x(j) contains the right side of equation (j)
c               as modified during elimination, j=1,...,last, while
c               for j .gt. last, x(j) is already a component of the
c               solution vector.
c       x(1),...,x(ncol) contains, on output, the components of the
c               solution corresponding to the present block.
c
c*********************************************************************
c
      integer  last,  i,j,k,km1,lm1,lp1
      double precision w(nrow,ncol),x(ncol), t
c
      lp1 = last + 1
      if ( lp1 .gt. ncol )                          go to 30
      do 20 j = lp1, ncol
         t = - x(j)
         if ( t .eq. 0.d0 )                         go to 20
         do 10 i = 1, last
   10       x(i) = x(i) + w(i,j) * t
   20 continue
   30 if ( last .eq. 1 )                            go to 60
      lm1 = last - 1
      do 50 kb = 1, lm1
        km1 = last - kb
        k = km1 + 1
        x(k) = x(k)/w(k,k)
        t = - x(k)
        if ( t .eq. 0.d0 )                          go to 50
        do 40 i = 1, km1
   40     x(i) = x(i) + w(i,k) * t
   50 continue
   60 x(1) = x(1)/w(1,1)
      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc timestamp prints out the current ymdhms date as a timestamp.
c
c  licensing:
c
c    this code is distributed under the gnu lgpl license.
c
c  modified:
c
c    12 january 2007
c
c  author:
c
c    john burkardt
c
c  parameters:
c
c    none
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'january  ', 'february ', 'march    ', 'april    ',
     &  'may      ', 'june     ', 'july     ', 'august   ',
     &  'september', 'october  ', 'november ', 'december ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'am'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'noon'
        else
          ampm = 'pm'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'pm'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'midnight'
          else
            ampm = 'am'
          end if
        end if
      end if

      write ( *,
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
