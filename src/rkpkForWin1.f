   
      DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DTEMP
      INTEGER I,INCX,M,MP1,N,NINCX
C
      DASUM = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DTEMP = DTEMP + DABS(DX(I))
   10 CONTINUE
      DASUM = DTEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I + 1)) + DABS(DX(I + 2))
     *  + DABS(DX(I + 3)) + DABS(DX(I + 4)) + DABS(DX(I + 5))
   50 CONTINUE
   60 DASUM = DTEMP
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
      subroutine dbimdr (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y,
&      tol1, tol2, init, prec1, maxiter1, prec2, maxiter2, theta, nlaht,
&      score, varht, c, d, eta, wk, swk, qwk, ywk, u, w, info)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxiter1, 
&     maxiter2, info
      double precision s(lds,*), q(ldqr,ldqc,*), y(2,*), tol1, tol2, 
&     prec1, prec2, theta(*), nlaht, score, varht, c(*), d(*), wk(*), 
&     eta(*), swk(lds,*), qwk(ldqr,ldqc,*), ywk(*), u(*), w(*)
c      character*2 vmu
      integer vmu
      double precision mse, tmp, dasum, mtol
      integer i, j
      info = 0
      mtol = 1.d0
23000 if(.not.( 1.d0 + mtol .gt. 1.d0 ))goto 23001
      mtol = mtol / 2.d0
      goto 23000
23001 continue
      if(.not.( mtol .lt. tol1 ))goto 23002
      mtol = tol1
23002 continue
23004 continue
      maxiter2 = maxiter2 - 1
      j=1
23007 if(.not.(j.le.nobs))goto 23009
      if(.not.(eta(j) .gt. 700.d0))goto 23010
      tmp = 1.d0
      goto 23011
23010 continue
      tmp = dexp (eta(j)) / (1.d0 + dexp (eta(j)))
23011 continue
      u(j) = y(1,j) * tmp - y(2,j)
      w(j) = y(1,j) * tmp * (1 - tmp)
      if(.not.(w(j) .le. mtol))goto 23012
      info = -7
      goto 23009
23012 continue
      i=1
23014 if(.not.(i.le.nnull))goto 23016
      swk(j,i) = s(j,i) * dsqrt (w(j))
      i=i+1
      goto 23014
23016 continue
      ywk(j) = dsqrt (w(j)) * (eta(j) - u(j) / w(j))
      j=j+1
      goto 23007
23009 continue
      if(.not.(info .eq. -7))goto 23017
      goto 23006
23017 continue
      call dcopy (ldqr*ldqc*nq, q, 1, qwk, 1)
      i=1
23019 if(.not.(i.le.nq))goto 23021
      j=1
23022 if(.not.(j.le.ldqr))goto 23024
      call dscal (ldqr-j+1, dsqrt (w(j)), qwk(j,j,i), 1)
      call dscal (j, dsqrt (w(j)), qwk(j,1,i), ldqr)
      j=j+1
      goto 23022
23024 continue
      i=i+1
      goto 23019
23021 continue
c      if(.not.(vmu .eq. 'u~'))
      if(.not.(vmu .eq. 3))goto 23025
      varht = 0.d0
      vmu=2
      j=1
23027 if(.not.(j.le.nobs))goto 23029
      varht = varht + u(j)**2 / w(j)
      j=j+1
      goto 23027
23029 continue
      varht = varht / dble (nobs)
23025 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dmudr (vmu, swk, lds, nobs, nnull, qwk, ldqr, ldqc, nq, ywk, 
&     tol2, init, prec1, maxiter1, theta, nlaht, score, varht, c, d, wk,
&      info)
      init = 1
      mse = 0.d0
      j=1
23030 if(.not.(j.le.nobs))goto 23032
      tmp = eta(j)
      eta(j) = (u(j) - 10.d0 ** nlaht * c(j)) / dsqrt (w(j))
      c(j) = c(j) * dsqrt (w(j))
      mse = mse + w(j) * ((eta(j)-tmp)/(1.d0+eta(j))) ** 2
      j=j+1
      goto 23030
23032 continue
      mse = dsqrt (mse / dasum (nobs, w, 1))
      if(.not.( info .ne. 0 ))goto 23033
      goto 23006
23033 continue
      if(.not.( mse .lt. prec2 ))goto 23035
      goto 23006
23035 continue
      if(.not.( maxiter2 .lt. 1 ))goto 23037
      info = -6
      goto 23006
23037 continue
23005 goto 23004
23006 continue
      return
      end
      subroutine dbisdr (vmu, s, lds, nobs, nnull, y, q, ldq, tol1, 
&     tol2, job, limnla, prec, maxiter, nlaht, score, varht, c, d, eta, 
&     qraux, jpvt, wk, swk, qwk, ywk, u, w, info)
c      character*2 vmu
      integer vmu
      integer lds, nobs, nnull, ldq, job, jpvt(*), info, maxiter
      double precision s(lds,*), y(2,*), q(ldq,*), tol1, tol2, limnla(2)
&     , nlaht, score(*), varht, c(*), d(*), qraux(*), wk(*), prec, eta(*
&     ), swk(lds,*), qwk(ldq,*), ywk(*), u(*), w(*)
      double precision mse, tmp, dasum, mtol
      integer i, j
      info = 0
      mtol = 1.d0
23000 if(.not.( 1.d0 + mtol .gt. 1.d0 ))goto 23001
      mtol = mtol / 2.d0
      goto 23000
23001 continue
      if(.not.( mtol .lt. tol1 ))goto 23002
      mtol = tol1
23002 continue
23004 continue
      maxiter = maxiter - 1
      j=1
23007 if(.not.(j.le.nobs))goto 23009
      if(.not.(eta(j) .gt. 700.d0))goto 23010
      tmp = 1.d0
      goto 23011
23010 continue
      tmp = dexp (eta(j)) / (1.d0 + dexp (eta(j)))
23011 continue
      u(j) = y(1,j) * tmp - y(2,j)
      w(j) = y(1,j) * tmp * (1 - tmp)
      if(.not.(w(j) .le. mtol))goto 23012
      info = -5
      goto 23009
23012 continue
      i=1
23014 if(.not.(i.le.nnull))goto 23016
      swk(j,i) = s(j,i) * dsqrt (w(j))
      i=i+1
      goto 23014
23016 continue
      ywk(j) = dsqrt (w(j)) * (eta(j) - u(j) / w(j))
      j=j+1
      goto 23007
23009 continue
      if(.not.(info .eq. -5))goto 23017
      goto 23006
23017 continue
      call dcopy (ldq*nobs, q, 1, qwk, 1)
      j=1
23019 if(.not.(j.le.nobs))goto 23021
      call dscal (nobs-j+1, dsqrt (w(j)), qwk(j,j), 1)
      call dscal (j, dsqrt (w(j)), qwk(j,1), nobs)
      j=j+1
      goto 23019
23021 continue
c      if(.not.(vmu .eq. 'u~'))
      if(.not.(vmu .eq. 3))goto 23022
      varht = 0.d0
      vmu=2
      j=1
23024 if(.not.(j.le.nobs))goto 23026
      varht = varht + u(j)**2 / w(j)
      j=j+1
      goto 23024
23026 continue
      varht = varht / dble (nobs)
23022 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dsidr (vmu, swk, lds, nobs, nnull, ywk, qwk, ldq, tol2, job, 
&     limnla, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
      mse = 0.d0
      j=1
23027 if(.not.(j.le.nobs))goto 23029
      tmp = eta(j)
      eta(j) = (u(j) - 10.d0 ** nlaht * c(j)) / dsqrt (w(j))
      c(j) = c(j) * dsqrt (w(j))
      mse = mse + w(j) * ((eta(j)-tmp)/(1.d0+dabs(eta(j)))) ** 2
      j=j+1
      goto 23027
23029 continue
      mse = dsqrt (mse / dasum (nobs, w, 1))
      if(.not.( info .ne. 0 ))goto 23030
      goto 23006
23030 continue
      if(.not.( mse .lt. prec ))goto 23032
      goto 23006
23032 continue
      if(.not.( maxiter .lt. 1 ))goto 23034
      info = -4
      goto 23006
23034 continue
23005 goto 23004
23006 continue
      return
      end
      subroutine dbmdr (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, 
&     tol1, tol2, init, prec1, maxiter1, prec2, maxiter2, theta, nlaht, 
&     score, varht, c, d, eta, wk, swk, qwk, ywk, u, w, info)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxiter1, 
&     maxiter2, info
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol1, tol2, 
&     prec1, prec2, theta(*), nlaht, score, varht, c(*), d(*), wk(*), 
&     eta(*), swk(lds,*), qwk(ldqr,ldqc,*), ywk(*), u(*), w(*)
c      character*2 vmu
      integer vmu
      double precision mse, tmp, dasum, mtol
      integer i, j
      info = 0
      mtol = 1.d0
23000 if(.not.( 1.d0 + mtol .gt. 1.d0 ))goto 23001
      mtol = mtol / 2.d0
      goto 23000
23001 continue
      if(.not.( mtol .lt. tol1 ))goto 23002
      mtol = tol1
23002 continue
23004 continue
      maxiter2 = maxiter2 - 1
      j=1
23007 if(.not.(j.le.nobs))goto 23009
      if(.not.(eta(j) .gt. 700.d0))goto 23010
      tmp = 1.d0
      goto 23011
23010 continue
      tmp = dexp (eta(j)) / (1.d0 + dexp (eta(j)))
23011 continue
      u(j) = tmp - y(j)
      w(j) = tmp * (1 - tmp)
      if(.not.(w(j) .le. mtol))goto 23012
      info = -7
      goto 23009
23012 continue
      i=1
23014 if(.not.(i.le.nnull))goto 23016
      swk(j,i) = s(j,i) * dsqrt (w(j))
      i=i+1
      goto 23014
23016 continue
      ywk(j) = dsqrt (w(j)) * (eta(j) - u(j) / w(j))
      j=j+1
      goto 23007
23009 continue
      if(.not.(info .eq. -7))goto 23017
      goto 23006
23017 continue
      call dcopy (ldqr*ldqc*nq, q, 1, qwk, 1)
      i=1
23019 if(.not.(i.le.nq))goto 23021
      j=1
23022 if(.not.(j.le.ldqr))goto 23024
      call dscal (ldqr-j+1, dsqrt (w(j)), qwk(j,j,i), 1)
      call dscal (j, dsqrt (w(j)), qwk(j,1,i), ldqr)
      j=j+1
      goto 23022
23024 continue
      i=i+1
      goto 23019
23021 continue
c      if(.not.(vmu .eq. 'u~'))
      if(.not.(vmu .eq. 3))goto 23025
      varht = 0.d0
      vmu=2
      j=1
23027 if(.not.(j.le.nobs))goto 23029
      varht = varht + u(j)**2 / w(j)
      j=j+1
      goto 23027
23029 continue
      varht = varht / dble (nobs)
23025 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dmudr (vmu, swk, lds, nobs, nnull, qwk, ldqr, ldqc, nq, ywk, 
&     tol2, init, prec1, maxiter1, theta, nlaht, score, varht, c, d, wk,
&      info)
      init = 1
      mse = 0.d0
      j=1
23030 if(.not.(j.le.nobs))goto 23032
      tmp = eta(j)
      eta(j) = (u(j) - 10.d0 ** nlaht * c(j)) / dsqrt (w(j))
      c(j) = c(j) * dsqrt (w(j))
      mse = mse + w(j) * ((eta(j)-tmp)/(1.d0+eta(j))) ** 2
      j=j+1
      goto 23030
23032 continue
      mse = dsqrt (mse / dasum (nobs, w, 1))
      if(.not.( info .ne. 0 ))goto 23033
      goto 23006
23033 continue
      if(.not.( mse .lt. prec2 ))goto 23035
      goto 23006
23035 continue
      if(.not.( maxiter2 .lt. 1 ))goto 23037
      info = -6
      goto 23006
23037 continue
23005 goto 23004
23006 continue
      return
      end
      subroutine dbsdr (vmu, s, lds, nobs, nnull, y, q, ldq, tol1, tol2,
&      job, limnla, prec, maxiter, nlaht, score, varht, c, d, eta, 
&     qraux, jpvt, wk, swk, qwk, ywk, u, w, info)
c      character*2 vmu
      integer vmu
      integer lds, nobs, nnull, ldq, job, jpvt(*), info, maxiter
      double precision s(lds,*), y(*), q(ldq,*), tol1, tol2, limnla(2), 
&     nlaht, score(*), varht(2), c(*), d(*), qraux(*), wk(*), prec, eta(
&     *), swk(lds,*), qwk(ldq,*), ywk(*), u(*), w(*)
      double precision mse, tmp, dasum, mtol
      integer i, j
      info = 0
      mtol = 1.d0
23000 if(.not.( 1.d0 + mtol .gt. 1.d0 ))goto 23001
      mtol = mtol / 2.d0
      goto 23000
23001 continue
      if(.not.( mtol .lt. tol1 ))goto 23002
      mtol = tol1
23002 continue
23004 continue
      maxiter = maxiter - 1
      j=1
23007 if(.not.(j.le.nobs))goto 23009
      if(.not.(eta(j) .gt. 700.d0))goto 23010
      tmp = 1.d0
      goto 23011
23010 continue
      tmp = dexp (eta(j)) / (1.d0 + dexp (eta(j)))
23011 continue
      u(j) = tmp - y(j)
      w(j) = tmp * (1 - tmp)
      if(.not.(w(j) .le. mtol))goto 23012
      info = -5
      goto 23009
23012 continue
      if(.not.( nnull .gt. 0 ))goto 23014
      i=1
23016 if(.not.(i.le.nnull))goto 23018
      swk(j,i) = s(j,i) * dsqrt (w(j))
      i=i+1
      goto 23016
23018 continue
23014 continue
      ywk(j) = dsqrt (w(j)) * (eta(j) - u(j) / w(j))
      j=j+1
      goto 23007
23009 continue
      if(.not.(info .eq. -5))goto 23019
      goto 23006
23019 continue
      call dcopy (ldq*nobs, q, 1, qwk, 1)
      j=1
23021 if(.not.(j.le.nobs))goto 23023
      call dscal (nobs-j+1, dsqrt (w(j)), qwk(j,j), 1)
      call dscal (j, dsqrt (w(j)), qwk(j,1), nobs)
      j=j+1
      goto 23021
23023 continue
c      if(.not.(vmu .eq. 'u~')) 
      if(.not.(vmu .eq. 3))goto 23024
      varht(1) = 0.d0
      vmu=2
      j=1
23026 if(.not.(j.le.nobs))goto 23028
      varht(1) = varht(1) + u(j)**2 / w(j)
      j=j+1
      goto 23026
23028 continue
      varht(1) = varht(1) / dble (nobs)
23024 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dsidr (vmu, swk, lds, nobs, nnull, ywk, qwk, ldq, tol2, job, 
&     limnla, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
      mse = 0.d0
      j=1
23029 if(.not.(j.le.nobs))goto 23031
      tmp = eta(j)
      eta(j) = (u(j) - 10.d0 ** nlaht * c(j)) / dsqrt (w(j))
      c(j) = c(j) * dsqrt (w(j))
      mse = mse + w(j) * ((eta(j)-tmp)/(1.d0+dabs(eta(j)))) ** 2
      j=j+1
      goto 23029
23031 continue
      mse = dsqrt (mse / dasum (nobs, w, 1))
      if(.not.( info .ne. 0 ))goto 23032
      goto 23006
23032 continue
      if(.not.( mse .lt. prec ))goto 23034
      goto 23006
23034 continue
      if(.not.( maxiter .lt. 1 ))goto 23036
      info = -4
      goto 23006
23036 continue
23005 goto 23004
23006 continue
      return
      end
      SUBROUTINE  DCOPY(N,DX,INCX,DY,INCY)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
      END
      subroutine ddeev (vmu, nobs, q, ldqr, ldqc, n, nq, u, ldu, uaux, 
&     t, x, theta, nlaht, score, varht, hes, ldh, gra, hwk1, hwk2, gwk1,
&      gwk2, kwk, ldk, work1, work2, work3, info)
c      character*1 vmu
      integer vmu
      integer nobs, ldqr, ldqc, n, nq, ldu, ldh, ldk, info
      double precision q(ldqr,ldqc,*), u(ldu,*), uaux(*), t(2,*), x(*), 
&     theta(*), nlaht, score, varht(2), hes(ldh,*), gra(*), hwk1(nq,*), 
&     hwk2(nq,*), gwk1(*), gwk2(*), kwk(ldk,ldk,*), work1(*), work2(*), 
&     work3(*)
      double precision trc, det, dum, ddot
      integer i, j, m
      info = 0
      call dset (nq, 0.d0, gra, 1)
      call dset (nq*nq, 0.d0, hes, 1)
c      if(.not.( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' ))
      if(.not.( vmu .ne. 0 .and. vmu .ne. 1 .and. vmu .ne. 2 ))
&     goto 23000
      info = -3
      return
23000 continue
      if(.not.( nobs .lt. n .or. ldqr .lt. n .or. ldqc .lt. n .or. nq 
&     .le. 0 .or. ldu .lt. n-1 .or. ldh .lt. nq .or. ldk .lt. n ))goto 2
&     3002
      info = -1
      return
23002 continue
      i=2
23004 if(.not.(i.le.nq))goto 23006
      if(.not.( theta(i) .le. -25.d0 ))goto 23007
      goto 23005
23007 continue
      j=1
23009 if(.not.(j.le.n))goto 23011
      call dcopy (n-j+1, q(j,j,i), 1, kwk(j,j,i), 1)
      call dscal (n-j+1, 10.d0 ** theta(i), kwk(j,j,i), 1)
      j=j+1
      goto 23009
23011 continue
      call dqrslm (u, ldu, n-1, n-2, uaux, kwk(2,2,i), n, 0, info, 
&     work1)
      call dqrsl (u, ldu, n-1, n-2, uaux, kwk(2,1,i), dum, kwk(2,1,i), 
&     dum, dum, dum, 01000, info)
23005 i=i+1
      goto 23004
23006 continue
      call dcopy (n, t(2,1), 2, kwk(1,1,1), n+1)
      call dcopy (n-1, t(1,2), 2, kwk(2,1,1), n+1)
      j=1
23012 if(.not.(j.lt.n-1))goto 23014
      call dset (n-j-1, 0.d0, kwk(j+2,j,1), 1)
      j=j+1
      goto 23012
23014 continue
      i=2
23015 if(.not.(i.le.nq))goto 23017
      if(.not.( theta(i) .le. -25.d0 ))goto 23018
      goto 23016
23018 continue
      j=1
23020 if(.not.(j.le.n))goto 23022
      call daxpy (n-j+1, -1.d0, kwk(j,j,i), 1, kwk(j,j,1), 1)
      j=j+1
      goto 23020
23022 continue
23016 i=i+1
      goto 23015
23017 continue
      i=1
23023 if(.not.(i.le.nq))goto 23025
      if(.not.( theta(i) .le. -25.d0 ))goto 23026
      goto 23024
23026 continue
      j=1
23028 if(.not.(j.lt.n))goto 23030
      call dcopy (n-j, kwk(j+1,j,i), 1, kwk(j,j+1,i), n)
      j=j+1
      goto 23028
23030 continue
23024 i=i+1
      goto 23023
23025 continue
      call dset (n, 10.d0 ** nlaht, work1, 1)
      call daxpy (n, 1.d0, work1, 1, t(2,1), 2)
      call dpbfa (t, 2, n, 1, info)
      if(.not.( info .ne. 0 ))goto 23031
      info = -2
      return
23031 continue
      i=1
23033 if(.not.(i.le.nq))goto 23035
      if(.not.( theta(i) .le. -25.d0 ))goto 23036
      goto 23034
23036 continue
      j=1
23038 if(.not.(j.le.n))goto 23040
      call dpbsl (t, 2, n, 1, kwk(1,j,i))
      j=j+1
      goto 23038
23040 continue
23034 i=i+1
      goto 23033
23035 continue
      call dcopy (n, x, 1, work1, 1)
      call dpbsl (t, 2, n, 1, work1)
c      if(.not.( vmu .ne. 'm' ))
      if(.not.( vmu .ne. 1 ))goto 23041
      call dcopy (n, work1, 1, work2, 1)
      call dscal (n, 2.d0, work2, 1)
      goto 23042
23041 continue
      call dcopy (n, x, 1, work2, 1)
23042 continue
      i=1
23043 if(.not.(i.le.nq))goto 23045
      if(.not.( theta(i) .le. -25.d0 ))goto 23046
      goto 23044
23046 continue
      call dgemv ('t', n, n, 1.d0, kwk(1,1,i), n, work2, 1, 0.d0, work3,
&      1)
      gwk1(i) = - ddot (n, work1, 1, work3, 1)
23044 i=i+1
      goto 23043
23045 continue
      i=1
23048 if(.not.(i.le.nq))goto 23050
      gwk2(i) = 0.d0
      if(.not.( theta(i) .le. -25.d0 ))goto 23051
      goto 23049
23051 continue
      j=1
23053 if(.not.(j.le.n))goto 23055
c      if(.not.( vmu .ne. 'm' ))
      if(.not.( vmu .ne. 1 ))goto 23056
      call dcopy (n, kwk(1,j,i), 1, work1, 1)
      call dpbsl (t, 2, n, 1, work1)
      gwk2(i) = gwk2(i) - work1(j)
      goto 23057
23056 continue
      gwk2(i) = gwk2(i) - kwk(j,j,i)
23057 continue
      j=j+1
      goto 23053
23055 continue
23049 i=i+1
      goto 23048
23050 continue
c      if(.not.( vmu .ne. 'm' ))
      if(.not.( vmu .ne. 1 ))goto 23058
      call dcopy (n, x, 1, work1, 1)
      call dpbsl (t, 2, n, 1, work1)
      i=1
23060 if(.not.(i.le.nq))goto 23062
      if(.not.( theta(i) .le. -25.d0 ))goto 23063
      goto 23061
23063 continue
      call dgemv ('n', n, n, 1.d0, kwk(1,1,i), n, work1, 1, 0.d0, work2,
&      1)
      j=1
23065 if(.not.(j.le.i))goto 23067
      if(.not.( theta(j) .le. -25.d0 ))goto 23068
      goto 23066
23068 continue
      call dgemv ('n', n, n, 1.d0, kwk(1,1,j), n, work1, 1, 0.d0, work3,
&      1)
      hwk1(i,j) = 2.d0 * ddot (n, work2, 1, work3, 1)
      call dgemv ('t', n, n, 1.d0, kwk(1,1,j), n, work1, 1, 0.d0, work3,
&      1)
      hwk1(i,j) = hwk1(i,j) + 2.d0 * ddot (n, work2, 1, work3, 1)
23066 j=j+1
      goto 23065
23067 continue
      call dgemv ('t', n, n, 1.d0, kwk(1,1,i), n, work1, 1, 0.d0, work2,
&      1)
      j=1
23070 if(.not.(j.le.i))goto 23072
      if(.not.( theta(j) .le. -25.d0 ))goto 23073
      goto 23071
23073 continue
      call dgemv ('n', n, n, 1.d0, kwk(1,1,j), n, work1, 1, 0.d0, work3,
&      1)
      hwk1(i,j) = hwk1(i,j) + 2.d0 * ddot (n, work2, 1, work3, 1)
23071 j=j+1
      goto 23070
23072 continue
23061 i=i+1
      goto 23060
23062 continue
      goto 23059
23058 continue
      call dcopy (n, x, 1, work1, 1)
      call dpbsl (t, 2, n, 1, work1)
      i=1
23075 if(.not.(i.le.nq))goto 23077
      if(.not.( theta(i) .le. -25.d0 ))goto 23078
      goto 23076
23078 continue
      call dgemv ('n', n, n, 1.d0, kwk(1,1,i), n, work1, 1, 0.d0, work2,
&      1)
      j=1
23080 if(.not.(j.le.i))goto 23082
      if(.not.( theta(j) .le. -25.d0 ))goto 23083
      goto 23081
23083 continue
      call dgemv ('t', n, n, 1.d0, kwk(1,1,j), n, x, 1, 0.d0, work3, 1)
      hwk1(i,j) = 2.d0 * ddot (n, work2, 1, work3, 1)
23081 j=j+1
      goto 23080
23082 continue
23076 i=i+1
      goto 23075
23077 continue
23059 continue
      i=1
23085 if(.not.(i.le.nq))goto 23087
      if(.not.( theta(i) .le. -25.d0 ))goto 23088
      goto 23086
23088 continue
      hwk1(i,i) = hwk1(i,i) + gwk1(i)
23086 i=i+1
      goto 23085
23087 continue
      i=1
23090 if(.not.(i.le.nq))goto 23092
      if(.not.( theta(i) .le. -25.d0 ))goto 23093
      goto 23091
23093 continue
      m=1
23095 if(.not.(m.le.i))goto 23097
      hwk2(i,m) = 0.d0
      if(.not.( theta(m) .le. -25.d0 ))goto 23098
      goto 23096
23098 continue
      j=1
23100 if(.not.(j.le.n))goto 23102
c      if(.not.( vmu .ne. 'm' ))
      if(.not.( vmu .ne. 1 ))goto 23103
      call dcopy (n, kwk(1,j,m), 1, work1, 1)
      call dpbsl (t, 2, n, 1, work1)
      hwk2(i,m) = hwk2(i,m) + 2.d0 * ddot (n, kwk(j,1,i), n, work1, 1)
      goto 23104
23103 continue
      hwk2(i,m) = hwk2(i,m) + ddot (n, kwk(j,1,i), n, kwk(1,j,m), 1)
23104 continue
      j=j+1
      goto 23100
23102 continue
23096 m=m+1
      goto 23095
23097 continue
23091 i=i+1
      goto 23090
23092 continue
      i=1
23105 if(.not.(i.le.nq))goto 23107
      if(.not.( theta(i) .le. -25.d0 ))goto 23108
      goto 23106
23108 continue
      hwk2(i,i) = hwk2(i,i) + gwk2(i)
23106 i=i+1
      goto 23105
23107 continue
c      if(.not.( vmu .eq. 'v' ))
      if(.not.( vmu .eq. 0 ))goto 23110
      trc = dfloat (nobs) * 10.d0 ** (-nlaht) * varht(1) / score
      i=1
23112 if(.not.(i.le.nq))goto 23114
      if(.not.( theta(i) .le. -25.d0 ))goto 23115
      goto 23113
23115 continue
      gra(i) = gwk1(i) / trc / trc - 2.d0 * score * gwk2(i) / trc / 
&     dfloat(nobs)
23113 i=i+1
      goto 23112
23114 continue
      call dscal (nq, dfloat (nobs), gra, 1)
23110 continue
c      if(.not.( vmu .eq. 'u' ))
      if(.not.( vmu .eq. 2 ))goto 23117
      dum = 10.d0 ** nlaht
      i=1
23119 if(.not.(i.le.nq))goto 23121
      if(.not.( theta(i) .le. -25.d0 ))goto 23122
      goto 23120
23122 continue
      gra(i) = dum * dum * gwk1(i) - 2.d0 * varht(1) * dum * gwk2(i)
23120 i=i+1
      goto 23119
23121 continue
      call dscal (nq, 1.d0/dfloat (n), gra, 1)
23117 continue
c      if(.not.( vmu .eq. 'm' ))
      if(.not.( vmu .eq. 1 ))goto 23124
      det = 10.d0 ** (-nlaht) * varht(2) / score
      i=1
23126 if(.not.(i.le.nq))goto 23128
      if(.not.( theta(i) .le. -25.d0 ))goto 23129
      goto 23127
23129 continue
      gra(i) = gwk1(i) / det - dfloat (nobs) / dfloat (n) * score * 
&     gwk2(i)
23127 i=i+1
      goto 23126
23128 continue
      call dscal (nq, 1.d0 / dfloat (nobs), gra, 1)
23124 continue
c      if(.not.( vmu .eq. 'v' ))
      if(.not.( vmu .eq. 0 ))goto 23131
      i=1
23133 if(.not.(i.le.nq))goto 23135
      if(.not.( theta(i) .le. -25.d0 ))goto 23136
      goto 23134
23136 continue
      j=1
23138 if(.not.(j.le.i))goto 23140
      if(.not.( theta(j) .le. -25.d0 ))goto 23141
      goto 23139
23141 continue
      hes(i,j) = hwk1(i,j) / trc / trc - 2.d0 * gwk1(i) * gwk2(j) / trc 
&     ** 3 - 2.d0 * gwk1(j) * gwk2(i) / trc ** 3 - 2.d0 * score * hwk2(
&     i,j) / trc / dfloat (nobs) + 6.d0 * score * gwk2(i) * gwk2(j) / 
&     trc / trc / dfloat (nobs)
23139 j=j+1
      goto 23138
23140 continue
      call dscal (i, dfloat (nobs), hes(i,1), ldh)
23134 i=i+1
      goto 23133
23135 continue
23131 continue
c      if(.not.( vmu .eq. 'u' ))
      if(.not.( vmu .eq. 2 ))goto 23143
      i=1
23145 if(.not.(i.le.nq))goto 23147
      if(.not.( theta(i) .le. -25.d0 ))goto 23148
      goto 23146
23148 continue
      j=1
23150 if(.not.(j.le.i))goto 23152
      if(.not.( theta(j) .le. -25.d0 ))goto 23153
      goto 23151
23153 continue
      hes(i,j) = dum * dum * hwk1(i,j) - 2.d0 * varht(1) * dum * hwk2(i,
&     j)
23151 j=j+1
      goto 23150
23152 continue
      call dscal (i, 1.d0/dfloat (n), hes(i,1), ldh)
23146 i=i+1
      goto 23145
23147 continue
23143 continue
c      if(.not.( vmu .eq. 'm' ))
      if(.not.( vmu .eq. 1 ))goto 23155
      i=1
23157 if(.not.(i.le.nq))goto 23159
      if(.not.( theta(i) .le. -25.d0 ))goto 23160
      goto 23158
23160 continue
      j=1
23162 if(.not.(j.le.i))goto 23164
      if(.not.( theta(j) .le. -25.d0 ))goto 23165
      goto 23163
23165 continue
      hes(i,j) = hwk1(i,j) / det - gwk1(i) * gwk2(j) / det / dfloat (n) 
&     - gwk1(j) * gwk2(i) / det / dfloat (n) - dfloat (nobs) / dfloat (
&     n) * score * hwk2(i,j) + dfloat (nobs) / dfloat (n) ** 2 * score *
&      gwk2(i) * gwk2(j)
23163 j=j+1
      goto 23162
23164 continue
      call dscal (i, 1.d0 / dfloat (nobs), hes(i,1), ldh)
23158 i=i+1
      goto 23157
23159 continue
23155 continue
      return
      end
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
      subroutine deval (vmu, q, ldq, n, z, nint, low, upp, nlaht, score,
&      varht, info, twk, work)
c      character*1 vmu
      integer vmu
      integer ldq, n, nint, info
      double precision q(ldq,*), z(*), low, upp, nlaht, score(*), varht(
&     2), twk(2,*), work(*)
      double precision tmp, minscr, mlo, varhtwk(2)
      integer j
      info = 0
      if(.not.( upp .lt. low ))goto 23000
      mlo = low
      low = upp
      upp = mlo
23000 continue
c      if(.not.( (vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u')
      if(.not.( (vmu .ne. 0 .and. vmu .ne. 1 .and. vmu .ne. 2)  
&     .or. nint .lt. 1 ))goto 23002
      info = -3
      return
23002 continue
      if(.not.( 1 .gt. n .or. n .gt. ldq ))goto 23004
      info = -1
      return
23004 continue
      j=1
23006 if(.not.(j.le.nint+1))goto 23008
      tmp = low + dfloat (j-1) * ( upp - low ) / dfloat (nint)
      call dset (n, 10.d0 ** (tmp), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**tmp
      call dtrev (vmu, twk, 2, n, z, score(j), varht, info, work)
      if(.not.( info .ne. 0 ))goto 23009
      info = -2
      return
23009 continue
      if(.not.( score(j) .le. minscr .or. j .eq. 1 ))goto 23011
      minscr = score(j)
      nlaht = tmp
      varhtwk(1) = varht(1)
      varhtwk(2) = varht(2)
23011 continue
      j=j+1
      goto 23006
23008 continue
      varht(1) = varhtwk(1)
      varht(2) = varhtwk(2)
      return
      end
      DOUBLE PRECISION FUNCTION devlpl(a,n,x)
C
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DEVLPL(A,N,X)
C              Double precision EVALuate a PoLynomial at X
C
C
C                              Function
C
C
C     returns
C          A(1) + A(2)*X + ... + A(N)*X**(N-1)
C
C
C                              Arguments
C
C
C     A --> Array of coefficients of the polynomial.
C                                        A is DOUBLE PRECISION(N)
C
C     N --> Length of A, also degree of polynomial - 1.
C                                        N is INTEGER
C
C     X --> Point at which the polynomial is to be evaluated.
C                                        X is DOUBLE PRECISION
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
      INTEGER n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION a(n)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION term
      INTEGER i
C     ..
C     .. Executable Statements ..
      term = a(n)
      DO 10,i = n - 1,1,-1
          term = a(i) + term*x
   10 CONTINUE
      devlpl = term
      RETURN

      END
************************************************************************
*
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
              JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
*
      subroutine dgmdr (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, 
&     tol1, tol2, init, prec1, maxiter1, prec2, maxiter2, theta, nlaht, 
&     score, varht, c, d, eta, wk, swk, qwk, ywk, u, w, info)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxiter1, 
&     maxiter2, info
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol1, tol2, 
&     prec1, prec2, theta(*), nlaht, score, varht, c(*), d(*), wk(*), 
&     eta(*), swk(lds,*), qwk(ldqr,ldqc,*), ywk(*), u(*), w(*)
c      character*2 vmu
      integer vmu
      double precision mse, tmp, dasum, mtol
      integer i, j
      info = 0
      mtol = 1.d0
23000 if(.not.( 1.d0 + mtol .gt. 1.d0 ))goto 23001
      mtol = mtol / 2.d0
      goto 23000
23001 continue
      if(.not.( mtol .lt. tol1 ))goto 23002
      mtol = tol1
23002 continue
23004 continue
      maxiter2 = maxiter2 - 1
      j=1
23007 if(.not.(j.le.nobs))goto 23009
      if(.not.(eta(j) .lt. -700.d0))goto 23010
      tmp = 1.d0
      goto 23011
23010 continue
      tmp = dexp (-eta(j))
23011 continue
      u(j) = 1.d0 - y(j) * tmp
      w(j) = y(j) * tmp
      if(.not.(w(j) .le. mtol))goto 23012
      info = -7
      goto 23009
23012 continue
      i=1
23014 if(.not.(i.le.nnull))goto 23016
      swk(j,i) = s(j,i) * dsqrt (w(j))
      i=i+1
      goto 23014
23016 continue
      ywk(j) = dsqrt (w(j)) * (eta(j) - u(j) / w(j))
      j=j+1
      goto 23007
23009 continue
      if(.not.(info .eq. -7))goto 23017
      goto 23006
23017 continue
      call dcopy (ldqr*ldqc*nq, q, 1, qwk, 1)
      i=1
23019 if(.not.(i.le.nq))goto 23021
      j=1
23022 if(.not.(j.le.ldqr))goto 23024
      call dscal (ldqr-j+1, dsqrt (w(j)), qwk(j,j,i), 1)
      call dscal (j, dsqrt (w(j)), qwk(j,1,i), ldqr)
      j=j+1
      goto 23022
23024 continue
      i=i+1
      goto 23019
23021 continue
c      if(.not.(vmu .eq. 'u~'))
      if(.not.(vmu .eq. 3))goto 23025
      varht = 0.d0
      vmu=2
      j=1
23027 if(.not.(j.le.nobs))goto 23029
      varht = varht + u(j)**2 / w(j)
      j=j+1
      goto 23027
23029 continue
      varht = varht / dble (nobs)
23025 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dmudr (vmu, swk, lds, nobs, nnull, qwk, ldqr, ldqc, nq, ywk, 
&     tol2, init, prec1, maxiter1, theta, nlaht, score, varht, c, d, wk,
&      info)
      init = 1
      mse = 0.d0
      j=1
23030 if(.not.(j.le.nobs))goto 23032
      tmp = eta(j)
      eta(j) = (u(j) - 10.d0 ** nlaht * c(j)) / dsqrt (w(j))
      c(j) = c(j) * dsqrt (w(j))
      mse = mse + w(j) * ((eta(j)-tmp)/(1.d0+eta(j))) ** 2
      j=j+1
      goto 23030
23032 continue
      mse = dsqrt (mse / dasum (nobs, w, 1))
      if(.not.( info .ne. 0 ))goto 23033
      goto 23006
23033 continue
      if(.not.( mse .lt. prec2 ))goto 23035
      goto 23006
23035 continue
      if(.not.( maxiter2 .lt. 1 ))goto 23037
      info = -6
      goto 23006
23037 continue
23005 goto 23004
23006 continue
      return
      end
      subroutine dgold (vmu, q, ldq, n, z, low, upp, nlaht, score, 
&     varht, info, twk, work)
c      character*1 vmu
      integer vmu
      integer ldq, n, info
      double precision q(ldq,*), z(*), low, upp, nlaht, score, varht(2),
&      twk(2,*), work(*)
      double precision ratio, mlo, mup, tmpl, tmpu
      ratio = ( dsqrt (5.d0) - 1.d0 ) / 2.d0
      info = 0
      if(.not.( upp .lt. low ))goto 23000
      mlo = low
      low = upp
      upp = mlo
23000 continue
c      if(.not.( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' ))
      if(.not.( vmu .ne. 0 .and. vmu .ne. 1 .and. vmu .ne. 2 ))
&     goto 23002
      info = -3
      return
23002 continue
      if(.not.( n .lt. 1 .or. n .gt. ldq ))goto 23004
      info = -1
      return
23004 continue
      mlo = upp - ratio * (upp - low)
      call dset (n, 10.d0 ** (mlo), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mlo
      call dtrev (vmu, twk, 2, n, z, tmpl, varht, info, work)
      if(.not.( info .ne. 0 ))goto 23006
      info = -2
      return
23006 continue
      mup = low + ratio * (upp - low)
      call dset (n, 10.d0 ** (mup), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mup
      call dtrev (vmu, twk, 2, n, z, tmpu, varht, info, work)
      if(.not.( info .ne. 0 ))goto 23008
      info = -2
      return
23008 continue
23010 continue
      if(.not.( mup - mlo .lt. 1.d-7 ))goto 23013
      goto 23012
23013 continue
      if(.not.( tmpl .lt. tmpu ))goto 23015
      upp = mup
      mup = mlo
      tmpu = tmpl
      mlo = upp - ratio * (upp - low)
      call dset (n, 10.d0 ** (mlo), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mlo
      call dtrev (vmu, twk, 2, n, z, tmpl, varht, info, work)
      if(.not.( info .ne. 0 ))goto 23017
      info = -2
      return
23017 continue
      goto 23016
23015 continue
      low = mlo
      mlo = mup
      tmpl = tmpu
      mup = low + ratio * (upp - low)
      call dset (n, 10.d0 ** (mup), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**mup
      call dtrev (vmu, twk, 2, n, z, tmpu, varht, info, work)
      if(.not.( info .ne. 0 ))goto 23019
      info = -2
      return
23019 continue
23016 continue
23011 goto 23010
23012 continue
      nlaht = ( mup + mlo ) / 2.d0
      call dset (n, 10.d0 ** (nlaht), twk(2,1), 2)
      call daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
      twk(1,1) = 10.d0**nlaht
      call dtrev (vmu, twk, 2, n, z, score, varht, info, work)
      if(.not.( info .ne. 0 ))goto 23021
      info = -2
      return
23021 continue
      return
      end
      subroutine dgsdr (vmu, s, lds, nobs, nnull, y, q, ldq, tol1, tol2,
&      job, limnla, prec, maxiter, nlaht, score, varht, c, d, eta, 
&     qraux, jpvt, wk, swk, qwk, ywk, u, w, info)
c      character*2 vmu
      integer vmu
      integer lds, nobs, nnull, ldq, job, jpvt(*), info, maxiter
      double precision s(lds,*), y(*), q(ldq,*), tol1, tol2, limnla(2), 
&     nlaht, score(*), varht, c(*), d(*), qraux(*), wk(*), prec, eta(*),
&      swk(lds,*), qwk(ldq,*), ywk(*), u(*), w(*)
      double precision mse, tmp, dasum, mtol
      integer i, j
      info = 0
      mtol = 1.d0
23000 if(.not.( 1.d0 + mtol .gt. 1.d0 ))goto 23001
      mtol = mtol / 2.d0
      goto 23000
23001 continue
      if(.not.( mtol .lt. tol1 ))goto 23002
      mtol = tol1
23002 continue
23004 continue
      maxiter = maxiter - 1
      j=1
23007 if(.not.(j.le.nobs))goto 23009
      if(.not.(eta(j) .lt. -700.d0))goto 23010
      tmp = 1.d0
      goto 23011
23010 continue
      tmp = dexp (-eta(j))
23011 continue
      u(j) = 1.d0 - y(j) * tmp
      w(j) = y(j) * tmp
      if(.not.(w(j) .le. mtol))goto 23012
      info = -5
      goto 23009
23012 continue
      i=1
23014 if(.not.(i.le.nnull))goto 23016
      swk(j,i) = s(j,i) * dsqrt (w(j))
      i=i+1
      goto 23014
23016 continue
      ywk(j) = dsqrt (w(j)) * (eta(j) - u(j) / w(j))
      j=j+1
      goto 23007
23009 continue
      if(.not.(info .eq. -5))goto 23017
      goto 23006
23017 continue
      call dcopy (ldq*nobs, q, 1, qwk, 1)
      j=1
23019 if(.not.(j.le.nobs))goto 23021
      call dscal (nobs-j+1, dsqrt (w(j)), qwk(j,j), 1)
      call dscal (j, dsqrt (w(j)), qwk(j,1), nobs)
      j=j+1
      goto 23019
23021 continue
c      if(.not.(vmu .eq. 'u~'))
      if(.not.(vmu .eq. 3))goto 23022
      varht = 0.d0
      vmu=2
      j=1
23024 if(.not.(j.le.nobs))goto 23026
      varht = varht + u(j)**2 / w(j)
      j=j+1
      goto 23024
23026 continue
      varht = varht / dble (nobs)
23022 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dsidr (vmu, swk, lds, nobs, nnull, ywk, qwk, ldq, tol2, job, 
&     limnla, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
      mse = 0.d0
      j=1
23027 if(.not.(j.le.nobs))goto 23029
      tmp = eta(j)
      eta(j) = (u(j) - 10.d0 ** nlaht * c(j)) / dsqrt (w(j))
      c(j) = c(j) * dsqrt (w(j))
      mse = mse + w(j) * ((eta(j)-tmp)/(1.d0+dabs(eta(j)))) ** 2
      j=j+1
      goto 23027
23029 continue
      mse = dsqrt (mse / dasum (nobs, w, 1))
      if(.not.( info .ne. 0 ))goto 23030
      goto 23006
23030 continue
      if(.not.( mse .lt. prec ))goto 23032
      goto 23006
23032 continue
      if(.not.( maxiter .lt. 1 ))goto 23034
      info = -4
      goto 23006
23034 continue
23005 goto 23004
23006 continue
      return
      end
      DOUBLE PRECISION FUNCTION dlanor(x)
C**********************************************************************
C
C      DOUBLE PRECISION FUNCTION DLANOR( X )
C           Double precision Logarith of the Asymptotic Normal
C
C
C                              Function
C
C
C      Computes the logarithm of the cumulative normal distribution
C      from abs( x ) to infinity for abs( x ) >= 5.
C
C
C                              Arguments
C
C
C      X --> Value at which cumulative normal to be evaluated
C                     DOUBLE PRECISION X
C
C
C                              Method
C
C
C      23 term expansion of formula 26.2.12 of Abramowitz and Stegun.
C      The relative error at X = 5 is about 0.5E-5.
C
C
C                              Note
C
C
C      ABS(X) must be >= 5 else there is an error stop.
C
C**********************************************************************
C     .. Parameters ..
      DOUBLE PRECISION dlsqpi
      PARAMETER (dlsqpi=0.91893853320467274177D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION approx,correc,xx,xx2
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION coef(12)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION devlpl,dln1px
      EXTERNAL devlpl,dln1px
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,log
C     ..
C     .. Data statements ..
      DATA coef/-1.0D0,3.0D0,-15.0D0,105.0D0,-945.0D0,10395.0D0,
     +     -135135.0D0,2027025.0D0,-34459425.0D0,654729075.0D0,
     +     -13749310575D0,316234143225.0D0/
C     ..
C     .. Executable Statements ..

      xx = abs(x)
      IF (xx.LT.5.0D0) STOP ' Argument too small in DLANOR'

      approx = -dlsqpi - 0.5*xx*xx - log(xx)

      xx2 = xx*xx
      correc = devlpl(coef,12,1.0D0/xx2)/xx2
      correc = dln1px(correc)

      dlanor = approx + correc

      RETURN

      END
      DOUBLE PRECISION FUNCTION dln1px(a)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DLN1PX(X)
C               Double precision LN(1+X)
C
C
C                              Function
C
C
C     Returns ln(1+x)
C     Note that the obvious code of
C               LOG(1.0+X)
C     won't work for small X because 1.0+X loses accuracy
C
C
C                              Arguments
C
C
C     X --> Value for which ln(1-x) is desired.
C                                        X is DOUBLE PRECISION
C
C
C                              Method
C
C
C     Renames ALNREL from:
C     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
C     Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
C     Trans. Math.  Softw. 18 (1993), 360-373.
C
C**********************************************************************
C-----------------------------------------------------------------------
C            EVALUATION OF THE FUNCTION LN(1 + A)
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION a
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION p1,p2,p3,q1,q2,q3,t,t2,w,x
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dlog
C     ..
C     .. Data statements ..
      DATA p1/-.129418923021993D+01/,p2/.405303492862024D+00/,
     +     p3/-.178874546012214D-01/
      DATA q1/-.162752256355323D+01/,q2/.747811014037616D+00/,
     +     q3/-.845104217945565D-01/
C     ..
C     .. Executable Statements ..
C--------------------------
      IF (abs(a).GT.0.375D0) GO TO 10
      t = a/ (a+2.0D0)
      t2 = t*t
      w = (((p3*t2+p2)*t2+p1)*t2+1.0D0)/ (((q3*t2+q2)*t2+q1)*t2+1.0D0)
      dln1px = 2.0D0*t*w
      RETURN
C
   10 x = 1.D0 + dble(a)
      dln1px = dlog(x)
      RETURN

      END
      subroutine dmcdc (a, lda, p, e, jpvt, info)
      integer lda, p, jpvt(*), info
      double precision a(lda,*), e(*)
      double precision beta, delta, theta, tmp, dasum, ddot
      integer i, j, jmax, jtmp, idamax
      info = 0
      if(.not.( lda .lt. p .or. p .lt. 1 ))goto 23000
      info = -1
      return
23000 continue
      tmp = 1.d0
23002 if(.not.( 1.d0 + tmp .gt. 1.d0 ))goto 23003
      tmp = tmp / 2.d0
      goto 23002
23003 continue
      jmax = idamax (p, a, lda+1)
      beta = dmax1 (2.d0 * tmp, dabs (a(jmax,jmax)))
      tmp = dsqrt (dfloat (p*p-1))
      if(.not.( tmp .lt. 1.d0 ))goto 23004
      tmp = 1.d0
23004 continue
      j=2
23006 if(.not.(j.le.p))goto 23008
      jmax = idamax (j-1, a(1,j), 1)
      beta = dmax1 (beta, dabs (a(jmax,j)) / tmp)
      j=j+1
      goto 23006
23008 continue
      delta = dasum (p, a, lda+1) / dfloat (p) * 1.d-7
      delta = dmax1 (delta, 1.d-10)
      j=1
23009 if(.not.(j.le.p))goto 23011
      jpvt(j) = j
      j=j+1
      goto 23009
23011 continue
      j=1
23012 if(.not.(j.le.p))goto 23014
      jmax = idamax (p-j+1, a(j,j), lda+1) + j - 1
      if(.not.( jmax .ne. j ))goto 23015
      call dswap (j-1, a(1,j), 1, a(1,jmax), 1)
      call dswap (jmax-j-1, a(j,j+1), lda, a(j+1,jmax), 1)
      call dswap (p-jmax, a(j,jmax+1), lda, a(jmax,jmax+1), lda)
      tmp = a(j,j)
      a(j,j) = a(jmax,jmax)
      a(jmax,jmax) = tmp
      jtmp = jpvt(j)
      jpvt(j) = jpvt(jmax)
      jpvt(jmax) = jtmp
23015 continue
      i=1
23017 if(.not.(i.lt.j))goto 23019
      a(i,j) = a(i,j) / a(i,i)
      i=i+1
      goto 23017
23019 continue
      i=j+1
23020 if(.not.(i.le.p))goto 23022
      a(j,i) = a(j,i) - ddot (j-1, a(1,j), 1, a(1,i), 1)
      i=i+1
      goto 23020
23022 continue
      if(.not.( j .eq. p ))goto 23023
      theta = 0.d0
      goto 23024
23023 continue
      jmax = idamax (p-j, a(j,j+1), lda) + j
      theta = dabs (a(j,jmax))
23024 continue
      tmp = dmax1 (delta, dabs (a(j,j)), theta ** 2 / beta)
      e(j) = tmp - a(j,j)
      a(j,j) = tmp
      i=j+1
23025 if(.not.(i.le.p))goto 23027
      a(i,i) = a(i,i) - a(j,i) ** 2 / a(j,j)
      i=i+1
      goto 23025
23027 continue
      j=j+1
      goto 23012
23014 continue
      j=1
23028 if(.not.(j.le.p))goto 23030
      a(j,j) = dsqrt (a(j,j))
      call dscal (p-j, a(j,j), a(j,j+1), lda)
      j=j+1
      goto 23028
23030 continue
      return
      end
      subroutine dmudr (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, 
&     tol, init, prec, maxite, theta, nlaht, score, varht, c, d, wk, 
&     info)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxite, info
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol, prec, theta(
&     *), nlaht, score, varht(2), c(*), d(*), wk(*)
c      character*1 vmu
      integer vmu
      integer n, n0
      integer iqraux, itraux, itwk, iqwk, iywk, ithewk, ihes, igra, 
&     ihwk1, ihwk2, igwk1, igwk2, ikwk, iwork1, iwork2, ijpvt, ipvtwk
      n = nobs
      n0 = nnull
      iqraux = 1
      itraux = iqraux + n0
      itwk = itraux + (n-n0-2)
      iqwk = itwk + 2 * (n-n0)
      iywk = iqwk + n * n
      ithewk = iywk + n
      ihes = ithewk + nq
      igra = ihes + nq * nq
      ihwk1 = igra + nq
      ihwk2 = ihwk1 + nq * nq
      igwk1 = ihwk2 + nq * nq
      igwk2 = igwk1 + nq
      ikwk = igwk2 + nq
      iwork1 = ikwk + (n-n0) * (n-n0) * nq
      iwork2 = iwork1 + n
      ijpvt = iwork2 + n
      ipvtwk = ijpvt + n0
      call dmudr1 (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, tol, 
&     init, prec, maxite, theta, nlaht, score, varht, c, d, wk(iqraux), 
&     wk(ijpvt), wk(itwk), wk(itraux), wk(iqwk), wk(iywk), wk(ithewk), 
&     wk(ihes), wk(igra), wk(ihwk1), wk(ihwk2), wk(igwk1), wk(igwk2), 
&     wk(ipvtwk), wk(ikwk), wk(iwork1), wk(iwork2), info)
      return
      end
      subroutine dmudrnew (vmu, s, lds, nobs, nnull, q,q1, q2, ldqr, 
&     ldqc,nq, y, tol, init, prec, maxite, theta, nlaht, score, varht, 
&     c, d, wk, ok, info)
      integer lds, nobs, nnull, ldqr, ldqc, init, maxite, info, nq, ok
      double precision s(lds,*),q(ldqr,ldqc,*), y(*), tol, prec, theta(*
&     ), nlaht, score, varht, c(*), d(*), wk(*), q1(ldqr,*), q2(ldqr,*)
c      character * 1 vmu
      integer vmu
      integer i, j
      i=1
23000 if(.not.(i.le.ldqr))goto 23002
      j=1
23003 if(.not.(j.le.ldqc))goto 23005
      q(i,j,1)=q1(i,j)
      q(i,j,2)=q2(i,j)
      j=j+1
      goto 23003
23005 continue
      i=i+1
      goto 23000
23002 continue
      call dmudr (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, tol, 
&     init, prec, maxite, theta, nlaht, score, varht, c, d, wk, info)
      if(.not.(d(1).lt. 0.001))goto 23006
      ok=1
23006 continue
      return
      end
      DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX)
      INTEGER          NEXT
      DOUBLE PRECISION   DX(1), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C
      IF(N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
      SUBROUTINE DPBFA(ABD,LDA,N,M,INFO)
      INTEGER LDA,N,M,INFO
      DOUBLE PRECISION ABD(LDA,1)
C
C     DPBFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
C     MATRIX STORED IN BAND FORM.
C
C     DPBFA IS USUALLY CALLED BY DPBCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C
C     ON ENTRY
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER
C                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE
C                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE
C                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. M + 1 .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        M       INTEGER
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. M .LT. N .
C
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND
C                FORM, SO THAT  A = TRANS(R)*R .
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT
C                     POSITIVE DEFINITE.
C
C     BAND STORAGE
C
C           IF  A  IS A SYMMETRIC POSITIVE DEFINITE BAND MATRIX,
C           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT.
C
C                   M = (BAND WIDTH ABOVE DIAGONAL)
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-M)
C                      DO 10 I = I1, J
C                         K = I-J+M+1
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DDOT
C     FORTRAN MAX0,DSQRT
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,T
      DOUBLE PRECISION S
      INTEGER IK,J,JK,K,MU
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
         DO 30 J = 1, N
            INFO = J
            S = 0.0D0
            IK = M + 1
            JK = MAX0(J-M,1)
            MU = MAX0(M+2-J,1)
            IF (M .LT. MU) GO TO 20
            DO 10 K = MU, M
               T = ABD(K,J) - DDOT(K-MU,ABD(IK,JK),1,ABD(MU,J),1)
               T = T/ABD(M+1,JK)
               ABD(K,J) = T
               S = S + T*T
               IK = IK - 1
               JK = JK + 1
   10       CONTINUE
   20       CONTINUE
            S = ABD(M+1,J) - S
C     ......EXIT
            IF (S .LE. 0.0D0) GO TO 40
            ABD(M+1,J) = DSQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
      SUBROUTINE DPBSL(ABD,LDA,N,M,B)
      INTEGER LDA,N,M
      DOUBLE PRECISION ABD(LDA,1),B(1)
C
C     DPBSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
C     BAND SYSTEM  A*X = B
C     USING THE FACTORS COMPUTED BY DPBCO OR DPBFA.
C
C     ON ENTRY
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DPBCO OR DPBFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        M       INTEGER
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
C        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
C        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
C        CORRECTLY AND  INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DPBCO(ABD,LDA,N,RCOND,Z,INFO)
C           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL DPBSL(ABD,LDA,N,C(1,J))
C        10 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C     FORTRAN MIN0
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,LA,LB,LM
C
C     SOLVE TRANS(R)*Y = B
C
      DO 10 K = 1, N
         LM = MIN0(K-1,M)
         LA = M + 1 - LM
         LB = K - LM
         T = DDOT(LM,ABD(LA,K),1,B(LB),1)
         B(K) = (B(K) - T)/ABD(M+1,K)
   10 CONTINUE
C
C     SOLVE R*X = Y
C
      DO 20 KB = 1, N
         K = N + 1 - KB
         LM = MIN0(K-1,M)
         LA = M + 1 - LM
         LB = K - LM
         B(K) = B(K)/ABD(M+1,K)
         T = -B(K)
         CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   20 CONTINUE
      RETURN
      END
      subroutine dpmdr (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, 
&     tol1, tol2, init, prec1, maxiter1, prec2, maxiter2, theta, nlaht, 
&     score, varht, c, d, eta, wk, swk, qwk, ywk, u, w, info)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxiter1, 
&     maxiter2, info
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol1, tol2, 
&     prec1, prec2, theta(*), nlaht, score, varht, c(*), d(*), wk(*), 
&     eta(*), swk(lds,*), qwk(ldqr,ldqc,*), ywk(*), u(*), w(*)
c      character*2 vmu
      integer vmu
      double precision mse, tmp, dasum, mtol
      integer i, j
      info = 0
      mtol = 1.d0
23000 if(.not.( 1.d0 + mtol .gt. 1.d0 ))goto 23001
      mtol = mtol / 2.d0
      goto 23000
23001 continue
      if(.not.( mtol .lt. tol1 ))goto 23002
      mtol = tol1
23002 continue
23004 continue
      maxiter2 = maxiter2 - 1
      j=1
23007 if(.not.(j.le.nobs))goto 23009
      if(.not.(eta(j) .gt. 700.d0))goto 23010
      w(j) = dexp(700.d0)
      goto 23011
23010 continue
      w(j) = dexp (eta(j))
23011 continue
      u(j) = w(j) - y(j)
      if(.not.(w(j) .le. mtol))goto 23012
      info = -7
      goto 23009
23012 continue
      i=1
23014 if(.not.(i.le.nnull))goto 23016
      swk(j,i) = s(j,i) * dsqrt (w(j))
      i=i+1
      goto 23014
23016 continue
      ywk(j) = dsqrt (w(j)) * (eta(j) - u(j) / w(j))
      j=j+1
      goto 23007
23009 continue
      if(.not.(info .eq. -7))goto 23017
      goto 23006
23017 continue
      call dcopy (ldqr*ldqc*nq, q, 1, qwk, 1)
      i=1
23019 if(.not.(i.le.nq))goto 23021
      j=1
23022 if(.not.(j.le.ldqr))goto 23024
      call dscal (ldqr-j+1, dsqrt (w(j)), qwk(j,j,i), 1)
      call dscal (j, dsqrt (w(j)), qwk(j,1,i), ldqr)
      j=j+1
      goto 23022
23024 continue
      i=i+1
      goto 23019
23021 continue
c      if(.not.(vmu .eq. 'u~'))
      if(.not.(vmu .eq. 3))goto 23025
      varht = 0.d0
      vmu=2
      j=1
23027 if(.not.(j.le.nobs))goto 23029
      varht = varht + u(j)**2 / w(j)
      j=j+1
      goto 23027
23029 continue
      varht = varht / dble (nobs)
23025 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dmudr (vmu, swk, lds, nobs, nnull, qwk, ldqr, ldqc, nq, ywk, 
&     tol2, init, prec1, maxiter1, theta, nlaht, score, varht, c, d, wk,
&      info)
      init = 1
      mse = 0.d0
      j=1
23030 if(.not.(j.le.nobs))goto 23032
      tmp = eta(j)
      eta(j) = (u(j) - 10.d0 ** nlaht * c(j)) / dsqrt (w(j))
      c(j) = c(j) * dsqrt (w(j))
      mse = mse + w(j) * ((eta(j)-tmp)/(1.d0+eta(j))) ** 2
      j=j+1
      goto 23030
23032 continue
      mse = dsqrt (mse / dasum (nobs, w, 1))
      if(.not.( info .ne. 0 ))goto 23033
      goto 23006
23033 continue
      if(.not.( mse .lt. prec2 ))goto 23035
      goto 23006
23035 continue
      if(.not.( maxiter2 .lt. 1 ))goto 23037
      info = -6
      goto 23006
23037 continue
23005 goto 23004
23006 continue
      return
      end
      SUBROUTINE DPOFA(A,LDA,N,INFO)
      INTEGER LDA,N,INFO
      DOUBLE PRECISION A(LDA,1)
C
C     DPOFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
C     MATRIX.
C
C     DPOFA IS USUALLY CALLED BY DPOCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DPOCO) = (1 + 18/N)*(TIME FOR DPOFA) .
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE
C                DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R
C                WHERE  TRANS(R)  IS THE TRANSPOSE.
C                THE STRICT LOWER TRIANGLE IS UNALTERED.
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DDOT
C     FORTRAN DSQRT
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,T
      DOUBLE PRECISION S
      INTEGER J,JM1,K
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
         DO 30 J = 1, N
            INFO = J
            S = 0.0D0
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 20
            DO 10 K = 1, JM1
               T = A(K,J) - DDOT(K-1,A(1,K),1,A(1,J),1)
               T = T/A(K,K)
               A(K,J) = T
               S = S + T*T
   10       CONTINUE
   20       CONTINUE
            S = A(J,J) - S
C     ......EXIT
            IF (S .LE. 0.0D0) GO TO 40
            A(J,J) = DSQRT(S)
   30    CONTINUE
         INFO = 0
   40 CONTINUE
      RETURN
      END
      SUBROUTINE DPOSL(A,LDA,N,B)
      INTEGER LDA,N
      DOUBLE PRECISION A(LDA,1),B(1)
C
C     DPOSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
C     SYSTEM A * X = B
C     USING THE FACTORS COMPUTED BY DPOCO OR DPOFA.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DPOCO OR DPOFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
C        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
C        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
C        CORRECTLY AND  INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DPOCO(A,LDA,N,RCOND,Z,INFO)
C           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
C           DO 10 J = 1, P
C              CALL DPOSL(A,LDA,N,C(1,J))
C        10 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB
C
C     SOLVE TRANS(R)*Y = B
C
      DO 10 K = 1, N
         T = DDOT(K-1,A(1,K),1,B(1),1)
         B(K) = (B(K) - T)/A(K,K)
   10 CONTINUE
C
C     SOLVE R*X = Y
C
      DO 20 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   20 CONTINUE
      RETURN
      END
      subroutine dprmut (x,npar,jpvt,job)
      integer npar,jpvt(npar),job
      double precision x(npar)
c
c Purpose: permute the elements of the array x according to the index 
c	vector jpvt (either forward or backward permutation).
c
c On Entry:
c   x(npar)		array to be permuted
c   npar		size of x (and jpvt)
c   jpvt		indices of the permutation
c   job			indicator of forward or backward permutation
c			if job = 0 forward permutation  
c				x(jpvt(i)) moved to x(i)
c			if job is nonzero backward permutation 
c				x(i) moved to x(jpvt(i))
c On Exit:
c   x(npar)		array with permuted entries
c
c   Written:	Yin Ling	U. of Maryland, August,1978
c
c $Header: dprmut.f,v 2.1 86/04/08 14:05:53 lindstrom Exp $
c
      integer i,j,k
      double precision t
c
      if (npar .le. 1) then
         return
      endif
      do 10 j = 1,npar
         jpvt(j) = -jpvt(j)
   10 continue
      if (job .eq. 0) then
c		forward permutation
         do 30 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 30
            endif
            j = i
            jpvt(j) = -jpvt(j)
            k = jpvt(j)
c           while
   20       if (jpvt(k) .lt. 0) then
               t = x(j)
               x(j) = x(k)
               x(k) = t
               jpvt(k) = -jpvt(k)
               j = k
               k = jpvt(k)
               goto 20
c           endwhile
            endif
   30    continue
      endif
      if (job .ne. 0 ) then
c			backward permutation
         do 50 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 50
            endif
            jpvt(i) = -jpvt(i)
            j = jpvt(i)
c           while
   40       if (j .ne. i) then
               t = x(i)
               x(i) = x(j)
               x(j) = t
               jpvt(j) = -jpvt(j)
               j = jpvt(j)
               goto 40
c           endwhile
            endif
   50    continue
      endif
      return
      end
      subroutine dpsdr (vmu, s, lds, nobs, nnull, y, q, ldq, tol1, tol2,
&      job, limnla, prec, maxiter, nlaht, score, varht, c, d, eta, 
&     qraux, jpvt, wk, swk, qwk, ywk, u, w, info)
c      character*2 vmu
      integer vmu
      integer lds, nobs, nnull, ldq, job, jpvt(*), info, maxiter
      double precision s(lds,*), y(*), q(ldq,*), tol1, tol2, limnla(2), 
&     nlaht, score(*), varht(2), c(*), d(*), qraux(*), wk(*), prec, eta(
&     *), swk(lds,*), qwk(ldq,*), ywk(*), u(*), w(*)
      double precision mse, tmp, dasum, mtol
      integer i, j
      info = 0
      mtol = 1.d0
23000 if(.not.( 1.d0 + mtol .gt. 1.d0 ))goto 23001
      mtol = mtol / 2.d0
      goto 23000
23001 continue
      if(.not.( mtol .lt. tol1 ))goto 23002
      mtol = tol1
23002 continue
23004 continue
      maxiter = maxiter - 1
      j=1
23007 if(.not.(j.le.nobs))goto 23009
      if(.not.(eta(j) .gt. 700.d0))goto 23010
      w(j) = dexp (700.d0)
      goto 23011
23010 continue
      w(j) = dexp (eta(j))
23011 continue
      u(j) = w(j) - y(j)
      if(.not.(w(j) .le. 0.1d0))goto 23012
      w(j) = 0.1d0
23012 continue
      if(.not.(w(j) .le. mtol))goto 23014
      info = -5
      goto 23009
23014 continue
      if(.not.( nnull .gt. 0 ))goto 23016
      i=1
23018 if(.not.(i.le.nnull))goto 23020
      swk(j,i) = s(j,i) * dsqrt (w(j))
      i=i+1
      goto 23018
23020 continue
23016 continue
      ywk(j) = dsqrt (w(j)) * (eta(j) - u(j) / w(j))
      j=j+1
      goto 23007
23009 continue
      if(.not.(info .eq. -5))goto 23021
      goto 23006
23021 continue
      call dcopy (ldq*nobs, q, 1, qwk, 1)
      j=1
23023 if(.not.(j.le.nobs))goto 23025
      call dscal (nobs-j+1, dsqrt (w(j)), qwk(j,j), 1)
      call dscal (j, dsqrt (w(j)), qwk(j,1), nobs)
      j=j+1
      goto 23023
23025 continue
c      if(.not.(vmu .eq. 'u~'))
      if(.not.(vmu .eq. 3))goto 23026
      varht(1) = 0.d0
      vmu=2
      j=1
23028 if(.not.(j.le.nobs))goto 23030
      varht(1) = varht(1) + u(j)**2 / w(j)
      j=j+1
      goto 23028
23030 continue
      varht(1) = varht(1) / dble (nobs)
23026 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dsidr (vmu, swk, lds, nobs, nnull, ywk, qwk, ldq, tol2, job, 
&     limnla, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
      mse = 0.d0
      j=1
23031 if(.not.(j.le.nobs))goto 23033
      tmp = eta(j)
      eta(j) = (u(j) - 10.d0 ** nlaht * c(j)) / dsqrt (w(j))
      c(j) = c(j) * dsqrt (w(j))
      mse = mse + w(j) * ((eta(j)-tmp)/(1.d0+dabs(eta(j)))) ** 2
      j=j+1
      goto 23031
23033 continue
      mse = dsqrt (mse / dasum (nobs, w, 1))
      if(.not.( info .ne. 0 ))goto 23034
      goto 23006
23034 continue
      if(.not.( mse .lt. prec ))goto 23036
      goto 23006
23036 continue
      if(.not.( maxiter .lt. 1 ))goto 23038
      info = -4
      goto 23006
23038 continue
23005 goto 23004
23006 continue
      return
      end
      SUBROUTINE DQRDC(X,LDX,N,P,QRAUX,JPVT,WORK,JOB)
      INTEGER LDX,N,P,JOB
      INTEGER JPVT(1)
      DOUBLE PRECISION X(LDX,1),QRAUX(1),WORK(1)
C
C     DQRDC USES HOUSEHOLDER TRANSFORMATIONS TO COMPUTE THE QR
C     FACTORIZATION OF AN N BY P MATRIX X.  COLUMN PIVOTING
C     BASED ON THE 2-NORMS OF THE REDUCED COLUMNS MAY BE
C     PERFORMED AT THE USERS OPTION.
C
C     ON ENTRY
C
C        X       DOUBLE PRECISION(LDX,P), WHERE LDX .GE. N.
C                X CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO BE
C                COMPUTED.
C
C        LDX     INTEGER.
C                LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C        N       INTEGER.
C                N IS THE NUMBER OF ROWS OF THE MATRIX X.
C
C        P       INTEGER.
C                P IS THE NUMBER OF COLUMNS OF THE MATRIX X.
C
C        JPVT    INTEGER(P).
C                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION
C                OF THE PIVOT COLUMNS.  THE K-TH COLUMN X(K) OF X
C                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE
C                VALUE OF JPVT(K).
C
C                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL
C                                      COLUMN.
C
C                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE COLUMN.
C
C                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL COLUMN.
C
C                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL COLUMNS
C                ARE MOVED TO THE BEGINNING OF THE ARRAY X AND FINAL
C                COLUMNS TO THE END.  BOTH INITIAL AND FINAL COLUMNS
C                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY
C                FREE COLUMNS ARE MOVED.  AT THE K-TH STAGE OF THE
C                REDUCTION, IF X(K) IS OCCUPIED BY A FREE COLUMN
C                IT IS INTERCHANGED WITH THE FREE COLUMN OF LARGEST
C                REDUCED NORM.  JPVT IS NOT REFERENCED IF
C                JOB .EQ. 0.
C
C        WORK    DOUBLE PRECISION(P).
C                WORK IS A WORK ARRAY.  WORK IS NOT REFERENCED IF
C                JOB .EQ. 0.
C
C        JOB     INTEGER.
C                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.
C                IF JOB .EQ. 0, NO PIVOTING IS DONE.
C                IF JOB .NE. 0, PIVOTING IS DONE.
C
C     ON RETURN
C
C        X       X CONTAINS IN ITS UPPER TRIANGLE THE UPPER
C                TRIANGULAR MATRIX R OF THE QR FACTORIZATION.
C                BELOW ITS DIAGONAL X CONTAINS INFORMATION FROM
C                WHICH THE ORTHOGONAL PART OF THE DECOMPOSITION
C                CAN BE RECOVERED.  NOTE THAT IF PIVOTING HAS
C                BEEN REQUESTED, THE DECOMPOSITION IS NOT THAT
C                OF THE ORIGINAL MATRIX X BUT THAT OF X
C                WITH ITS COLUMNS PERMUTED AS DESCRIBED BY JPVT.
C
C        QRAUX   DOUBLE PRECISION(P).
C                QRAUX CONTAINS FURTHER INFORMATION REQUIRED TO RECOVER
C                THE ORTHOGONAL PART OF THE DECOMPOSITION.
C
C        JPVT    JPVT(K) CONTAINS THE INDEX OF THE COLUMN OF THE
C                ORIGINAL MATRIX THAT HAS BEEN INTERCHANGED INTO
C                THE K-TH COLUMN, IF PIVOTING WAS REQUESTED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     DQRDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2
C     FORTRAN DABS,DMAX1,MIN0,DSQRT
C
C     INTERNAL VARIABLES
C
      INTEGER J,JP,L,LP1,LUP,MAXJ,PL,PU
      DOUBLE PRECISION MAXNRM,DNRM2,TT
      DOUBLE PRECISION DDOT,NRMXL,T
      LOGICAL NEGJ,SWAPJ
C
C
      PL = 1
      PU = 0
      IF (JOB .EQ. 0) GO TO 60
C
C        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS
C        ACCORDING TO JPVT.
C
         DO 20 J = 1, P
            SWAPJ = JPVT(J) .GT. 0
            NEGJ = JPVT(J) .LT. 0
            JPVT(J) = J
            IF (NEGJ) JPVT(J) = -J
            IF (.NOT.SWAPJ) GO TO 10
               IF (J .NE. PL) CALL DSWAP(N,X(1,PL),1,X(1,J),1)
               JPVT(J) = JPVT(PL)
               JPVT(PL) = J
               PL = PL + 1
   10       CONTINUE
   20    CONTINUE
         PU = P
         DO 50 JJ = 1, P
            J = P - JJ + 1
            IF (JPVT(J) .GE. 0) GO TO 40
               JPVT(J) = -JPVT(J)
               IF (J .EQ. PU) GO TO 30
                  CALL DSWAP(N,X(1,PU),1,X(1,J),1)
                  JP = JPVT(PU)
                  JPVT(PU) = JPVT(J)
                  JPVT(J) = JP
   30          CONTINUE
               PU = PU - 1
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
C
C     COMPUTE THE NORMS OF THE FREE COLUMNS.
C
      IF (PU .LT. PL) GO TO 80
      DO 70 J = PL, PU
         QRAUX(J) = DNRM2(N,X(1,J),1)
         WORK(J) = QRAUX(J)
   70 CONTINUE
   80 CONTINUE
C
C     PERFORM THE HOUSEHOLDER REDUCTION OF X.
C
      LUP = MIN0(N,P)
      DO 200 L = 1, LUP
         IF (L .LT. PL .OR. L .GE. PU) GO TO 120
C
C           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT
C           INTO THE PIVOT POSITION.
C
            MAXNRM = 0.0D0
            MAXJ = L
            DO 100 J = L, PU
               IF (QRAUX(J) .LE. MAXNRM) GO TO 90
                  MAXNRM = QRAUX(J)
                  MAXJ = J
   90          CONTINUE
  100       CONTINUE
            IF (MAXJ .EQ. L) GO TO 110
               CALL DSWAP(N,X(1,L),1,X(1,MAXJ),1)
               QRAUX(MAXJ) = QRAUX(L)
               WORK(MAXJ) = WORK(L)
               JP = JPVT(MAXJ)
               JPVT(MAXJ) = JPVT(L)
               JPVT(L) = JP
  110       CONTINUE
  120    CONTINUE
         QRAUX(L) = 0.0D0
         IF (L .EQ. N) GO TO 190
C
C           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L.
C
            NRMXL = DNRM2(N-L+1,X(L,L),1)
            IF (NRMXL .EQ. 0.0D0) GO TO 180
               IF (X(L,L) .NE. 0.0D0) NRMXL = DSIGN(NRMXL,X(L,L))
               CALL DSCAL(N-L+1,1.0D0/NRMXL,X(L,L),1)
               X(L,L) = 1.0D0 + X(L,L)
C
C              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
C              UPDATING THE NORMS.
C
               LP1 = L + 1
               IF (P .LT. LP1) GO TO 170
               DO 160 J = LP1, P
                  T = -DDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
                  CALL DAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
                  IF (J .LT. PL .OR. J .GT. PU) GO TO 150
                  IF (QRAUX(J) .EQ. 0.0D0) GO TO 150
                     TT = 1.0D0 - (DABS(X(L,J))/QRAUX(J))**2
                     TT = DMAX1(TT,0.0D0)
                     T = TT
                     TT = 1.0D0 + 0.05D0*TT*(QRAUX(J)/WORK(J))**2
                     IF (TT .EQ. 1.0D0) GO TO 130
                        QRAUX(J) = QRAUX(J)*DSQRT(T)
                     GO TO 140
  130                CONTINUE
                        QRAUX(J) = DNRM2(N-L,X(L+1,J),1)
                        WORK(J) = QRAUX(J)
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
  170          CONTINUE
C
C              SAVE THE TRANSFORMATION.
C
               QRAUX(L) = X(L,L)
               X(L,L) = -NRMXL
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
      RETURN
      END
      SUBROUTINE DQRSL(X,LDX,N,K,QRAUX,Y,QY,QTY,B,RSD,XB,JOB,INFO)
      INTEGER LDX,N,K,JOB,INFO
      DOUBLE PRECISION X(LDX,1),QRAUX(1),Y(1),QY(1),QTY(1),B(1),RSD(1),
     *                 XB(1)
C
C     DQRSL APPLIES THE OUTPUT OF DQRDC TO COMPUTE COORDINATE
C     TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS.
C     FOR K .LE. MIN(N,P), LET XK BE THE MATRIX
C
C            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))
C
C     FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL
C     N X P MATRIX X THAT WAS INPUT TO DQRDC (IF NO PIVOTING WAS
C     DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR
C     ORIGINAL ORDER).  DQRDC PRODUCES A FACTORED ORTHOGONAL MATRIX Q
C     AND AN UPPER TRIANGULAR MATRIX R SUCH THAT
C
C              XK = Q * (R)
C                       (0)
C
C     THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS
C     X AND QRAUX.
C
C     ON ENTRY
C
C        X      DOUBLE PRECISION(LDX,P).
C               X CONTAINS THE OUTPUT OF DQRDC.
C
C        LDX    INTEGER.
C               LDX IS THE LEADING DIMENSION OF THE ARRAY X.
C
C        N      INTEGER.
C               N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST
C               HAVE THE SAME VALUE AS N IN DQRDC.
C
C        K      INTEGER.
C               K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K
C               MUST NNOT BE GREATER THAN MIN(N,P), WHERE P IS THE
C               SAME AS IN THE CALLING SEQUENCE TO DQRDC.
C
C        QRAUX  DOUBLE PRECISION(P).
C               QRAUX CONTAINS THE AUXILIARY OUTPUT FROM DQRDC.
C
C        Y      DOUBLE PRECISION(N)
C               Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED
C               BY DQRSL.
C
C        JOB    INTEGER.
C               JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS
C               THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING
C               MEANING.
C
C                    IF A.NE.0, COMPUTE QY.
C                    IF B,C,D, OR E .NE. 0, COMPUTE QTY.
C                    IF C.NE.0, COMPUTE B.
C                    IF D.NE.0, COMPUTE RSD.
C                    IF E.NE.0, COMPUTE XB.
C
C               NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB
C               AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR
C               WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING
C               SEQUENCE.
C
C     ON RETURN
C
C        QY     DOUBLE PRECISION(N).
C               QY CONNTAINS Q*Y, IF ITS COMPUTATION HAS BEEN
C               REQUESTED.
C
C        QTY    DOUBLE PRECISION(N).
C               QTY CONTAINS TRANS(Q)*Y, IF ITS COMPUTATION HAS
C               BEEN REQUESTED.  HERE TRANS(Q) IS THE
C               TRANSPOSE OF THE MATRIX Q.
C
C        B      DOUBLE PRECISION(K)
C               B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM
C
C                    MINIMIZE NORM2(Y - XK*B),
C
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT
C               IF PIVOTING WAS REQUESTED IN DQRDC, THE J-TH
C               COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J)
C               OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO DQRDC.)
C
C        RSD    DOUBLE PRECISION(N).
C               RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B,
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS
C               ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE
C               ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK.
C
C        XB     DOUBLE PRECISION(N).
C               XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B,
C               IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO
C               THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE
C               OF X.
C
C        INFO   INTEGER.
C               INFO IS ZERO UNLESS THE COMPUTATION OF B HAS
C               BEEN REQUESTED AND R IS EXACTLY SINGULAR.  IN
C               THIS CASE, INFO IS THE INDEX OF THE FIRST ZERO
C               DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED.
C
C     THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED
C     IF THEIR COMPUTATION IS NOT REQUESTED AND IN THIS CASE
C     CAN BE REPLACED BY DUMMY VARIABLES IN THE CALLING PROGRAM.
C     TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME
C     ARRAY FOR DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A
C     FREQUENTLY OCCURING EXAMPLE IS WHEN ONE WISHES TO COMPUTE
C     ANY OF B, RSD, OR XB AND DOES NOT NEED Y OR QTY.  IN THIS
C     CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD, OR XB, WHILE
C     PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE
C     COMPUTED.  THUS THE CALLING SEQUENCE
C
C          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
C
C     WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD
C     OVERWRITING Y.  MORE GENERALLY, EACH ITEM IN THE FOLLOWING
C     LIST CONTAINS GROUPS OF PERMISSIBLE IDENTIFICATIONS FOR
C     A SINGLE CALLINNG SEQUENCE.
C
C          1. (Y,QTY,B) (RSD) (XB) (QY)
C
C          2. (Y,QTY,RSD) (B) (XB) (QY)
C
C          3. (Y,QTY,XB) (B) (RSD) (QY)
C
C          4. (Y,QY) (QTY,B) (RSD) (XB)
C
C          5. (Y,QY) (QTY,RSD) (B) (XB)
C
C          6. (Y,QY) (QTY,XB) (B) (RSD)
C
C     IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO
C     THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     DQRSL USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
C
C     BLAS DAXPY,DCOPY,DDOT
C     FORTRAN DABS,MIN0,MOD
C
C     INTERNAL VARIABLES
C
      INTEGER I,J,JJ,JU,KP1
      DOUBLE PRECISION DDOT,T,TEMP
      LOGICAL CB,CQY,CQTY,CR,CXB
C
C
C     SET INFO FLAG.
C
      INFO = 0
C
C     DETERMINE WHAT IS TO BE COMPUTED.
C
      CQY = JOB/10000 .NE. 0
      CQTY = MOD(JOB,10000) .NE. 0
      CB = MOD(JOB,1000)/100 .NE. 0
      CR = MOD(JOB,100)/10 .NE. 0
      CXB = MOD(JOB,10) .NE. 0
      JU = MIN0(K,N-1)
C
C     SPECIAL ACTION WHEN N=1.
C
      IF (JU .NE. 0) GO TO 40
         IF (CQY) QY(1) = Y(1)
         IF (CQTY) QTY(1) = Y(1)
         IF (CXB) XB(1) = Y(1)
         IF (.NOT.CB) GO TO 30
            IF (X(1,1) .NE. 0.0D0) GO TO 10
               INFO = 1
            GO TO 20
   10       CONTINUE
               B(1) = Y(1)/X(1,1)
   20       CONTINUE
   30    CONTINUE
         IF (CR) RSD(1) = 0.0D0
      GO TO 250
   40 CONTINUE
C
C        SET UP TO COMPUTE QY OR QTY.
C
         IF (CQY) CALL DCOPY(N,Y,1,QY,1)
         IF (CQTY) CALL DCOPY(N,Y,1,QTY,1)
         IF (.NOT.CQY) GO TO 70
C
C           COMPUTE QY.
C
            DO 60 JJ = 1, JU
               J = JU - JJ + 1
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 50
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -DDOT(N-J+1,X(J,J),1,QY(J),1)/X(J,J)
                  CALL DAXPY(N-J+1,T,X(J,J),1,QY(J),1)
                  X(J,J) = TEMP
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
         IF (.NOT.CQTY) GO TO 100
C
C           COMPUTE TRANS(Q)*Y.
C
            DO 90 J = 1, JU
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 80
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  T = -DDOT(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)
                  CALL DAXPY(N-J+1,T,X(J,J),1,QTY(J),1)
                  X(J,J) = TEMP
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        SET UP TO COMPUTE B, RSD, OR XB.
C
         IF (CB) CALL DCOPY(K,QTY,1,B,1)
         KP1 = K + 1
         IF (CXB) CALL DCOPY(K,QTY,1,XB,1)
         IF (CR .AND. K .LT. N) CALL DCOPY(N-K,QTY(KP1),1,RSD(KP1),1)
         IF (.NOT.CXB .OR. KP1 .GT. N) GO TO 120
            DO 110 I = KP1, N
               XB(I) = 0.0D0
  110       CONTINUE
  120    CONTINUE
         IF (.NOT.CR) GO TO 140
            DO 130 I = 1, K
               RSD(I) = 0.0D0
  130       CONTINUE
  140    CONTINUE
         IF (.NOT.CB) GO TO 190
C
C           COMPUTE B.
C
            DO 170 JJ = 1, K
               J = K - JJ + 1
               IF (X(J,J) .NE. 0.0D0) GO TO 150
                  INFO = J
C           ......EXIT
                  GO TO 180
  150          CONTINUE
               B(J) = B(J)/X(J,J)
               IF (J .EQ. 1) GO TO 160
                  T = -B(J)
                  CALL DAXPY(J-1,T,X(1,J),1,B,1)
  160          CONTINUE
  170       CONTINUE
  180       CONTINUE
  190    CONTINUE
         IF (.NOT.CR .AND. .NOT.CXB) GO TO 240
C
C           COMPUTE RSD OR XB AS REQUIRED.
C
            DO 230 JJ = 1, JU
               J = JU - JJ + 1
               IF (QRAUX(J) .EQ. 0.0D0) GO TO 220
                  TEMP = X(J,J)
                  X(J,J) = QRAUX(J)
                  IF (.NOT.CR) GO TO 200
                     T = -DDOT(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)
                     CALL DAXPY(N-J+1,T,X(J,J),1,RSD(J),1)
  200             CONTINUE
                  IF (.NOT.CXB) GO TO 210
                     T = -DDOT(N-J+1,X(J,J),1,XB(J),1)/X(J,J)
                     CALL DAXPY(N-J+1,T,X(J,J),1,XB(J),1)
  210             CONTINUE
                  X(J,J) = TEMP
  220          CONTINUE
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
      RETURN
      END
      subroutine dqrslm (x, ldx, n, k, qraux, a, lda, job, info, work)
      integer ldx, n, k, lda, job, info
      double precision x(ldx,*), qraux(*), a(lda,*), work(*)
      double precision tmp, alph, ddot
      integer i, j, step
      info = 0
      if(.not.( lda .lt. n .or. n .lt. k .or. k .lt. 1 ))goto 23000
      info = -1
      return
23000 continue
      if(.not.( job .ne. 0 .and. job .ne. 1 ))goto 23002
      info = 1
      return
23002 continue
      if(.not.( job .eq. 0 ))goto 23004
      j = 1
      step = 1
      goto 23005
23004 continue
      j = k
      step = -1
23005 continue
23006 if(.not.( j .ge. 1 .and. j .le. k ))goto 23007
      if(.not.( qraux(j) .eq. 0.0d0 ))goto 23008
      j = j + step
      goto 23006
23008 continue
      tmp = x(j,j)
      x(j,j) = qraux(j)
      i=1
23010 if(.not.(i.lt.j))goto 23012
      alph = - ddot (n-j+1, x(j,j), 1, a(j,i), 1) / x(j,j)
      call daxpy (n-j+1, alph, x(j,j), 1, a(j,i), 1)
      i=i+1
      goto 23010
23012 continue
      alph = 1.d0 / x(j,j)
      call dsymv ('l', n-j+1, alph, a(j,j), lda, x(j,j), 1, 0.d0, work(
     *j), 1)
      alph = - ddot (n-j+1, work(j), 1, x(j,j), 1) / 2.d0 / x(j,j)
      call daxpy (n-j+1, alph, x(j,j), 1, work(j), 1)
      call dsyr2 ('l', n-j+1, -1.d0, x(j,j), 1, work(j), 1, a(j,j), lda)
      x(j,j) = tmp
      j = j + step
      goto 23006
23007 continue
      return
      end
      SUBROUTINE  DSCAL(N,DA,DX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DA,DX(1)
      INTEGER I,INCX,M,MP1,N,NINCX
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
      subroutine  dset(n,da,dx,incx)
      integer n,incx
      double precision da,dx(*)
c
c Purpose : set vector dx to constant da. Unrolled loops are used for 
c	increment equal to one.
c
c On Entry:
c   n			length of dx
c   da			any constant
c   incx		increment for dx
c
c On Exit:
c   dx(n)		vector with all n entries set to da
c
c $Header: dset.f,v 2.1 86/04/08 14:06:25 lindstrom Exp $
c
      integer i,m,mp1,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da
        dx(i + 1) = da
        dx(i + 2) = da
        dx(i + 3) = da
        dx(i + 4) = da
   50 continue
      return
      end
      subroutine dsms (s, lds, nobs, nnull, jpvt, q, ldq, nlaht, sms, 
     *ldsms, wk, info)
      integer lds, nobs, nnull, jpvt(*), ldq, ldsms, info
      double precision s(lds,*), q(ldq,*), nlaht, sms(ldsms,*), wk(2,*)
      double precision dum, ddot
      integer i, j, n, n0
      info = 0
      if(.not.( nnull .lt. 1 .or. nnull .ge. nobs .or. nobs .gt. lds 
     *.or. nobs .gt. ldq .or. ldsms .lt. nnull ))goto 23000
      info = -1
      return
23000 continue
      n0 = nnull
      n = nobs - nnull
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, wk, 1)
      j=1
23002 if(.not.(j.le.n0))goto 23004
      call dcopy (n, q(n0+1,j), 1, q(j,n0+1), ldq)
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, wk, q(n0+2,j), dum, q(n0+
     *2,j), dum, dum, dum, 01000, info)
      j=j+1
      goto 23002
23004 continue
      call dset (n, 10.d0 ** nlaht, wk(2,1), 2)
      call daxpy (n, 1.d0, q(n0+1,n0+1), ldq+1, wk(2,1), 2)
      call dcopy (n-1, q(n0+1,n0+2), ldq+1, wk(1,2), 2)
      call dpbfa (wk, 2, n, 1, info)
      if(.not.( info .ne. 0 ))goto 23005
      info = -2
      return
23005 continue
      j=1
23007 if(.not.(j.le.n0))goto 23009
      call dpbsl (wk, 2, n, 1, q(n0+1,j))
      j=j+1
      goto 23007
23009 continue
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, wk, 1)
      j=1
23010 if(.not.(j.le.n0))goto 23012
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, wk, q(n0+2,j), q(n0+2,j),
     * dum, dum, dum, dum, 10000, info)
      j=j+1
      goto 23010
23012 continue
      i=1
23013 if(.not.(i.le.n0))goto 23015
      j=1
23016 if(.not.(j.lt.i))goto 23018
      sms(i,j) = sms(j,i)
      j=j+1
      goto 23016
23018 continue
      j=i
23019 if(.not.(j.le.n0))goto 23021
      sms(i,j) = q(j,i) - ddot (n, q(n0+1,j), 1, q(i,n0+1), ldq)
      j=j+1
      goto 23019
23021 continue
      sms(i,i) = sms(i,i) + 10.d0**nlaht
      i=i+1
      goto 23013
23015 continue
      j=1
23022 if(.not.(j.le.n0))goto 23024
      call dtrsl (s, lds, n0, sms(1,j), 01, info)
      j=j+1
      goto 23022
23024 continue
      i=1
23025 if(.not.(i.le.n0))goto 23027
      call dcopy (n0, sms(i,1), ldsms, wk, 1)
      call dtrsl (s, lds, n0, wk, 01, info)
      call dprmut (wk, n0, jpvt, 1)
      call dcopy (n0, wk, 1, sms(i,1), ldsms)
      i=i+1
      goto 23025
23027 continue
      j=1
23028 if(.not.(j.le.n0))goto 23030
      call dprmut (sms(1,j), n0, jpvt, 1)
      j=j+1
      goto 23028
23030 continue
      j=1
23031 if(.not.(j.le.n0))goto 23033
      call dcopy (n, q(j,n0+1), ldq, q(n0+1,j), 1)
      j=j+1
      goto 23031
23033 continue
      return
      end
c      this is a subroutine to sort a arrary, using heapsrot method
c      see "numerical recipes"
      
       subroutine dsort(n,ra)
       integer n
       double precision ra(n)

       l=n/2+1
       ir=n
 10    continue
           if (l .gt. 1) then
               l=l-1
	       rra=ra(l)
	   else
	       rra=ra(ir)
	       ra(ir)=ra(1)
	       ir=ir-1
	       if (ir .eq. 1) then
		   ra(1)=rra
		   return
	       endif
            endif	 		
	    i=l
	    j=l+l
 20         if (j .le. ir) then
                if (j .lt. ir) then
	            if (ra(j) .lt. ra(j+1)) j=j+1
                endif
	        if (rra .lt. ra(j)) then
		    ra(i)=ra(j)
		    i=j
		    j=j+j
		else
		    j=ir+1
		endif
 	    go to 20
	    endif
	    ra(i)=rra
       go to 10
       end

      subroutine dstup (s, lds, nobs, nnull, qraux, jpvt, y, q, ldqr, 
     *ldqc, nq, info, work)
      integer lds, nobs, nnull, jpvt(*), ldqr, ldqc, nq, info
      double precision s(lds,*), y(*), qraux(*), q(ldqr,ldqc,*), work(*)
      double precision dum
      integer j
      info = 0
      if(.not.( nobs .lt. 1 .or. nobs .gt. lds .or. nobs .gt. ldqr .or. 
     *nobs .gt. ldqc ))goto 23000
      info = -1
      return
23000 continue
      j=1
23002 if(.not.(j.le.nnull))goto 23004
      jpvt(j) = 0
      j=j+1
      goto 23002
23004 continue
      call dqrdc (s, lds, nobs, nnull, qraux, jpvt, work, 1)
      call dqrsl (s, lds, nobs, nnull, qraux, y, dum, y, work, dum, dum,
     * 01100, info)
      if(.not.( info .ne. 0 ))goto 23005
      return
23005 continue
      j=1
23007 if(.not.(j.le.nq))goto 23009
      call dqrslm (s, lds, nobs, nnull, qraux, q(1,1,j), ldqr, 0, info, 
     *work)
      j=j+1
      goto 23007
23009 continue
      return
      end
      SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
   50 CONTINUE
      RETURN
      END
************************************************************************
*
      SUBROUTINE DSYMV ( UPLO, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYMV  performs the matrix-vector  operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 5
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 7
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set up the start points in  X  and  Y.
*
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( N - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, N
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, N
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, N
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, N
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  y  when A is stored in upper triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               TEMP1 = ALPHA*X( J )
               TEMP2 = ZERO
               DO 50, I = 1, J - 1
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   50          CONTINUE
               Y( J ) = Y( J ) + TEMP1*A( J, J ) + ALPHA*TEMP2
   60       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 80, J = 1, N
               TEMP1 = ALPHA*X( JX )
               TEMP2 = ZERO
               IX    = KX
               IY    = KY
               DO 70, I = 1, J - 1
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
                  IX      = IX      + INCX
                  IY      = IY      + INCY
   70          CONTINUE
               Y( JY ) = Y( JY ) + TEMP1*A( J, J ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y  when A is stored in lower triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 100, J = 1, N
               TEMP1  = ALPHA*X( J )
               TEMP2  = ZERO
               Y( J ) = Y( J )       + TEMP1*A( J, J )
               DO 90, I = J + 1, N
                  Y( I ) = Y( I ) + TEMP1*A( I, J )
                  TEMP2  = TEMP2  + A( I, J )*X( I )
   90          CONTINUE
               Y( J ) = Y( J ) + ALPHA*TEMP2
  100       CONTINUE
         ELSE
            JX = KX
            JY = KY
            DO 120, J = 1, N
               TEMP1   = ALPHA*X( JX )
               TEMP2   = ZERO
               Y( JY ) = Y( JY )       + TEMP1*A( J, J )
               IX      = JX
               IY      = JY
               DO 110, I = J + 1, N
                  IX      = IX      + INCX
                  IY      = IY      + INCY
                  Y( IY ) = Y( IY ) + TEMP1*A( I, J )
                  TEMP2   = TEMP2   + A( I, J )*X( IX )
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP2
               JX      = JX      + INCX
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DSYMV .
*
      END
*
************************************************************************
*
      SUBROUTINE DSYR2 ( UPLO, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, N
      CHARACTER*1        UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYR2  performs the symmetric rank 2 operation
*
*     A := alpha*x*y' + alpha*y*x' + A,
*
*  where alpha is a scalar, x and y are n element vectors and A is an n
*  by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP1, TEMP2
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO, 'U' ).AND.
     $         .NOT.LSAME( UPLO, 'L' )      )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYR2 ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Set up the start points in X and Y if the increments are not both
*     unity.
*
      IF( ( INCX.NE.1 ).OR.( INCY.NE.1 ) )THEN
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( N - 1 )*INCX
         END IF
         IF( INCY.GT.0 )THEN
            KY = 1
         ELSE
            KY = 1 - ( N - 1 )*INCY
         END IF
         JX = KX
         JY = KY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF( LSAME( UPLO, 'U' ) )THEN
*
*        Form  A  when A is stored in the upper triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 20, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  DO 10, I = 1, J
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   10             CONTINUE
               END IF
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = KX
                  IY    = KY
                  DO 30, I = 1, J
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   30             CONTINUE
               END IF
               JX = JX + INCX
               JY = JY + INCY
   40       CONTINUE
         END IF
      ELSE
*
*        Form  A  when A is stored in the lower triangle.
*
         IF( ( INCX.EQ.1 ).AND.( INCY.EQ.1 ) )THEN
            DO 60, J = 1, N
               IF( ( X( J ).NE.ZERO ).OR.( Y( J ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( J )
                  TEMP2 = ALPHA*X( J )
                  DO 50, I = J, N
                     A( I, J ) = A( I, J ) + X( I )*TEMP1 + Y( I )*TEMP2
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( ( X( JX ).NE.ZERO ).OR.( Y( JY ).NE.ZERO ) )THEN
                  TEMP1 = ALPHA*Y( JY )
                  TEMP2 = ALPHA*X( JX )
                  IX    = JX
                  IY    = JY
                  DO 70, I = J, N
                     A( I, J ) = A( I, J ) + X( IX )*TEMP1
     $                                     + Y( IY )*TEMP2
                     IX        = IX        + INCX
                     IY        = IY        + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
               JY = JY + INCY
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DSYR2 .
*
      END
*
      subroutine dsytr (x, ldx, n, tol, info, work)
      integer ldx, n, info
      double precision x(ldx,*), tol, work(*)
      double precision nrmtot, nrmxj, alph, toltot, tolcum, toluni, dn, 
     *ddot
      integer j
      info = 0
      if(.not.( ldx .lt. n .or. n .le. 2 ))goto 23000
      info = -1
      return
23000 continue
      nrmtot = ddot (n, x, ldx+1, x, ldx+1)
      j=1 
23002 if(.not.(j.lt.n))goto 23004
      nrmtot = nrmtot + 2.d0 * ddot (n-j, x(j+1,j), 1, x(j+1,j), 1)
       j=j+1 
      goto 23002
23004 continue
      toltot = 1.d0
23005 if(.not.( 1.d0 + toltot .gt. 1.d0 ))goto 23006
      toltot = toltot / 2.d0
      goto 23005
23006 continue
      toltot = 4.d0 * toltot ** 2
      if(.not.( toltot .lt. tol ))goto 23007
      toltot = tol
23007 continue
      toltot = toltot * nrmtot
      dn = dfloat (n)
      toluni = toltot * 6.d0 / dn / ( dn - 1.d0 ) / ( 2.d0 * dn - 1.d0 )
      tolcum = 0.d0
      j=1 
23009 if(.not.(j.lt.n-1))goto 23011
      nrmtot = nrmtot - x(j,j) * x(j,j)
      nrmxj = ddot (n-j, x(j+1,j), 1, x(j+1,j), 1)
      dn = dfloat (n-j)
      tolcum = tolcum + toluni * dn * dn
      if(.not.( 2.d0 * nrmxj .le. tolcum ))goto 23012
      x(j,j+1) = 0.d0
      call dscal (n-j, 0.d0, x(j+1,j), 1)
      tolcum = tolcum - 2.d0 * nrmxj
      toltot = toltot - 2.d0 * nrmxj
      goto 23010
23012 continue
      if(.not.( x(j+1,j) .lt. 0.d0 ))goto 23014
      x(j,j+1) = dsqrt (nrmxj)
      goto 23015
23014 continue
      x(j,j+1) = - dsqrt (nrmxj)
23015 continue
      nrmtot = nrmtot - 2.d0 * nrmxj
      call dscal (n-j, -1.d0/x(j,j+1), x(j+1,j), 1)
      x(j+1,j) = 1.d0 + x(j+1,j)
      alph = 1.d0 / x(j+1,j)
      call dsymv ('l', n-j, alph, x(j+1,j+1), ldx, x(j+1,j), 1, 0.d0, 
     *work(j+1), 1)
      alph = - ddot (n-j, work(j+1), 1, x(j+1,j), 1) / 2.d0 / x(j+1,j)
      call daxpy (n-j, alph, x(j+1,j), 1, work(j+1), 1)
      call dsyr2 ('l', n-j, -1.d0, x(j+1,j), 1, work(j+1), 1, x(j+1,j+1)
     *, ldx)
23010  j=j+1 
      goto 23009
23011 continue
      x(n-1,n) = x(n,n-1)
      return
      end
      subroutine dtrev (vmu, t, ldt, n, z, score, varht, info, work)      
c      character*1 vmu
      integer vmu
      integer n, info
      double precision t(ldt,*), z(*), score, varht(2), work(*)
      double precision nume, deno, tmp, alph, la, dasum, ddot
      integer j
      info = 0
c      if(.not.( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' ))
      if(.not.( vmu .ne. 0 .and. vmu .ne. 1 .and. vmu .ne. 2 ))
&     goto 23000
      info = -3
      return
23000 continue
      la = t(1,1)
      alph = dfloat (n) / dasum (n, t(2,1), ldt)
      call dscal (n, alph, t(2,1), ldt)
      call dscal (n-1, alph, t(1,2), ldt)
      call dpbfa (t, ldt, n, 1, info)
      if(.not.( info .ne. 0 ))goto 23002
      return
23002 continue
      call dcopy (n, z, 1, work, 1)
      call dpbsl (t, ldt, n, 1, work)
c      if(.not.( vmu .eq. 'v' ))
      if(.not.( vmu .eq. 0 ))goto 23004
      tmp = 1.d0 / t(2,n) / t(2,n)
      deno = tmp
      j=n-1
23006 if(.not.(j.gt.0))goto 23008
      tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
      deno = deno + tmp
      j=j-1
      goto 23006
23008 continue
      nume = ddot (n, work, 1, work, 1) / dfloat (n)
      deno = deno / dfloat (n)
      varht(1) = alph * la * nume / deno
      score = nume / deno / deno
      deno = dlog (t(2,n))
      j=n-1
23009 if(.not.(j.gt.0))goto 23011
      deno = deno + dlog (t(2,j))
      j=j-1
      goto 23009
23011 continue
      nume = ddot (n, z, 1, work, 1) / dfloat (n)
      varht(2) = alph * la * nume
23004 continue
c      if(.not.( vmu .eq. 'm' ))
      if(.not.( vmu .eq. 1 ))goto 23012
      deno = dlog (t(2,n))
      j=n-1
23014 if(.not.(j.gt.0))goto 23016
      deno = deno + dlog (t(2,j))
      j=j-1
      goto 23014
23016 continue
      nume = ddot (n, z, 1, work, 1) / dfloat (n)
      varht(2) = alph * la * nume
      score = nume * dexp (2.d0 * deno / dfloat (n))
      tmp = 1.d0 / t(2,n) / t(2,n)
      deno = tmp
      j=n-1
23017 if(.not.(j.gt.0))goto 23019
      tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
      deno = deno + tmp
      j=j-1
      goto 23017
23019 continue
      nume = ddot (n, work, 1, work, 1) / dfloat (n)
      deno = deno / dfloat (n)
      varht(1) = alph * la * nume / deno
23012 continue
c      if(.not.( vmu .eq. 'u' ))
      if(.not.( vmu .eq. 2 ))goto 23020
      nume = ddot (n, work, 1, work, 1) / dfloat (n)
      tmp = 1.d0 / t(2,n) / t(2,n)
      deno = tmp
      j=n-1
23022 if(.not.(j.gt.0))goto 23024
      tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
      deno = deno + tmp
      j=j-1
      goto 23022
23024 continue
      deno = deno / dfloat (n)
      score = alph * alph * la * la * nume - 2.d0 * varht(1) * alph * 
&     la * deno
23020 continue
      return
      end
      SUBROUTINE DTRSL(T,LDT,N,B,JOB,INFO)
      INTEGER LDT,N,JOB,INFO
      DOUBLE PRECISION T(LDT,1),B(1)
C
C
C     DTRSL SOLVES SYSTEMS OF THE FORM
C
C                   T * X = B
C     OR
C                   TRANS(T) * X = B
C
C     WHERE T IS A TRIANGULAR MATRIX OF ORDER N. HERE TRANS(T)
C     DENOTES THE TRANSPOSE OF THE MATRIX T.
C
C     ON ENTRY
C
C         T         DOUBLE PRECISION(LDT,N)
C                   T CONTAINS THE MATRIX OF THE SYSTEM. THE ZERO
C                   ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
C                   THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
C                   USED TO STORE OTHER INFORMATION.
C
C         LDT       INTEGER
C                   LDT IS THE LEADING DIMENSION OF THE ARRAY T.
C
C         N         INTEGER
C                   N IS THE ORDER OF THE SYSTEM.
C
C         B         DOUBLE PRECISION(N).
C                   B CONTAINS THE RIGHT HAND SIDE OF THE SYSTEM.
C
C         JOB       INTEGER
C                   JOB SPECIFIES WHAT KIND OF SYSTEM IS TO BE SOLVED.
C                   IF JOB IS
C
C                        00   SOLVE T*X=B, T LOWER TRIANGULAR,
C                        01   SOLVE T*X=B, T UPPER TRIANGULAR,
C                        10   SOLVE TRANS(T)*X=B, T LOWER TRIANGULAR,
C                        11   SOLVE TRANS(T)*X=B, T UPPER TRIANGULAR.
C
C     ON RETURN
C
C         B         B CONTAINS THE SOLUTION, IF INFO .EQ. 0.
C                   OTHERWISE B IS UNALTERED.
C
C         INFO      INTEGER
C                   INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR.
C                   OTHERWISE INFO CONTAINS THE INDEX OF
C                   THE FIRST ZERO DIAGONAL ELEMENT OF T.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     G. W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DDOT
C     FORTRAN MOD
C
C     INTERNAL VARIABLES
C
      DOUBLE PRECISION DDOT,TEMP
      INTEGER CASE,J,JJ
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 150
C
C        CHECK FOR ZERO DIAGONAL ELEMENTS.
C
         DO 10 INFO = 1, N
C     ......EXIT
            IF (T(INFO,INFO) .EQ. 0.0D0) GO TO 150
   10    CONTINUE
         INFO = 0
C
C        DETERMINE THE TASK AND GO TO IT.
C
         CASE = 1
         IF (MOD(JOB,10) .NE. 0) CASE = 2
         IF (MOD(JOB,100)/10 .NE. 0) CASE = CASE + 2
         GO TO (20,50,80,110), CASE
C
C        SOLVE T*X=B FOR T LOWER TRIANGULAR
C
   20    CONTINUE
            B(1) = B(1)/T(1,1)
            IF (N .LT. 2) GO TO 40
            DO 30 J = 2, N
               TEMP = -B(J-1)
               CALL DAXPY(N-J+1,TEMP,T(J,J-1),1,B(J),1)
               B(J) = B(J)/T(J,J)
   30       CONTINUE
   40       CONTINUE
         GO TO 140
C
C        SOLVE T*X=B FOR T UPPER TRIANGULAR.
C
   50    CONTINUE
            B(N) = B(N)/T(N,N)
            IF (N .LT. 2) GO TO 70
            DO 60 JJ = 2, N
               J = N - JJ + 1
               TEMP = -B(J+1)
               CALL DAXPY(J,TEMP,T(1,J+1),1,B(1),1)
               B(J) = B(J)/T(J,J)
   60       CONTINUE
   70       CONTINUE
         GO TO 140
C
C        SOLVE TRANS(T)*X=B FOR T LOWER TRIANGULAR.
C
   80    CONTINUE
            B(N) = B(N)/T(N,N)
            IF (N .LT. 2) GO TO 100
            DO 90 JJ = 2, N
               J = N - JJ + 1
               B(J) = B(J) - DDOT(JJ-1,T(J+1,J),1,B(J+1),1)
               B(J) = B(J)/T(J,J)
   90       CONTINUE
  100       CONTINUE
         GO TO 140
C
C        SOLVE TRANS(T)*X=B FOR T UPPER TRIANGULAR.
C
  110    CONTINUE
            B(1) = B(1)/T(1,1)
            IF (N .LT. 2) GO TO 130
            DO 120 J = 2, N
               B(J) = B(J) - DDOT(J-1,T(1,J),1,B(1),1)
               B(J) = B(J)/T(J,J)
  120       CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION dumnor(x)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DUMNOR(X)
C
C
C                              Function
C
C
C     Computes the cumulative  of    the  normal   distribution,   i.e.,
C     the integral from -infinity to x of
C          (1/sqrt(2*pi)) exp(-u*u/2) du
C
C
C                              Method
C
C
C     The    rational function  approximation   from pages   90  - 92 of
C     Kennedy  and Gentle,  Statistical  Computing,  Marcel  Dekker,  NY
C     1980.
C
C
C                              Arguments
C
C
C     X --> Argument at which cumulative normal is evaluated
C                    DOUBLE PRECISION X
C
C**********************************************************************
C
C
C     PIM12 IS PI**(-1/2)
C     SQRT2 IS SQRT(2)
C
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION derf,derfc,pim12,sqrt2,z,z2,zm2
      LOGICAL qdirct
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION xden1(4),xden2(8),xden3(5),xnum1(4),xnum2(8),
     +                 xnum3(5)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION devlpl,dlanor
      EXTERNAL devlpl,dlanor
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp
C     ..
C     .. Data statements ..
      DATA xnum1/2.4266795523053175D2,2.1979261618294152D1,
     +     6.9963834886191355D0,-3.5609843701815385D-2/
      DATA xden1/2.1505887586986120D2,9.1164905404514901D1,
     +     1.5082797630407787D1,1.0000000000000000D0/
      DATA xnum2/3.004592610201616005D2,4.519189537118729422D2,
     +     3.393208167343436870D2,1.529892850469404039D2,
     +     4.316222722205673530D1,7.211758250883093659D0,
     +     5.641955174789739711D-1,-1.368648573827167067D-7/
      DATA xden2/3.004592609569832933D2,7.909509253278980272D2,
     +     9.313540948506096211D2,6.389802644656311665D2,
     +     2.775854447439876434D2,7.700015293522947295D1,
     +     1.278272731962942351D1,1.000000000000000000D0/
      DATA xnum3/-2.99610707703542174D-3,-4.94730910623250734D-2,
     +     -2.26956593539686930D-1,-2.78661308609647788D-1,
     +     -2.23192459734184686D-2/
      DATA xden3/1.06209230528467918D-2,1.91308926107829841D-1,
     +     1.05167510706793207D0,1.98733201817135256D0,
     +     1.00000000000000000D0/
      DATA pim12/0.5641895835477562869480795D0/
      DATA sqrt2/1.4142135623730950488D0/
C     ..
C     .. Executable Statements ..
      IF (.NOT. (abs(x).LT.1.0E-30)) GO TO 10
      dumnor = 0.5
      RETURN

      GO TO 50

   10 IF (.NOT. (x.LE.-38.0)) GO TO 20
      dumnor = 0.0
      RETURN

      GO TO 50

   20 IF (.NOT. (x.LE.-15.0)) GO TO 30
      dumnor = exp(dlanor(x))
      RETURN

      GO TO 50

   30 IF (.NOT. (x.GT.6.0)) GO TO 40
      dumnor = 1.0
      RETURN

      GO TO 50

   40 CONTINUE
   50 z = abs(x/sqrt2)
      z2 = z*z
      zm2 = 1.0D0/z2
      IF (z.LT.0.5D0) THEN
          derf = z*devlpl(xnum1,4,z2)/devlpl(xden1,4,z2)
          qdirct = .TRUE.

      ELSE IF (z.LT.4.0D0) THEN
          derfc = exp(-z2)*devlpl(xnum2,8,z)/devlpl(xden2,8,z)
          qdirct = .FALSE.

      ELSE
          derfc = (exp(-z2)/z)* (pim12+zm2*devlpl(xnum3,5,zm2)/
     +            devlpl(xden3,5,zm2))
          qdirct = .FALSE.
      END IF

      IF (.NOT. (x.GE.0.0)) GO TO 60
      IF (.NOT. (qdirct)) derf = 1.0D0 - derfc
      dumnor = (1.0D0+derf)/2.0D0
      GO TO 70

   60 IF (qdirct) derfc = 1.0D0 - derf
      dumnor = derfc/2.0D0
   70 RETURN

      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
C
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DMAX
      INTEGER I,INCX,IX,N
C
      IDAMAX = 0
      IF( N .LT. 1 ) RETURN
      IDAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5
         IDAMAX = I
         DMAX = DABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
         IF(DABS(DX(I)).LE.DMAX) GO TO 30
         IDAMAX = I
         DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END
      REAL FUNCTION invnor(p)
C
C**********************************************************************
C
C     REAL FUNCTION INVNOR(P)
C                    NORmal distribution INVerse
C
C
C                              Function
C
C
C     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
C     infinity to X of (1/SQRT(2*PI)) EXP(-U*U) dU is P
C
C
C                              Arguments
C
C
C     P --> The probability whose normal deviate is sought.
C                    P is REAL
C
C
C                              Method
C
C
C     The  rational   function   on  page 95    of Kennedy  and  Gentle,
C     Statistical Computing, Marcel Dekker, NY , 1980.
C
C
C                              Note
C
C
C     If P .lt. 1.0e-20 then INVNOR returns -10.0.
C     If P .ge. 1.0 then INVNOR returns 10.0.
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      REAL p
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION sign,y,z
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION xden(5),xnum(5)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION devlpl
      EXTERNAL devlpl
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,log,sqrt
C     ..
C     .. Data statements ..
      DATA xnum/-0.322232431088D0,-1.000000000000D0,-0.342242088547D0,
     +     -0.204231210245D-1,-0.453642210148D-4/
      DATA xden/0.993484626060D-1,0.588581570495D0,0.531103462366D0,
     +     0.103537752850D0,0.38560700634D-2/
C     ..
C     .. Executable Statements ..
      IF (.NOT. (p.LT.1.0E-20)) GO TO 10
      invnor = -10.0
      RETURN

   10 IF (.NOT. (p.GE.1.0)) GO TO 20
      invnor = 10.0
      RETURN

   20 IF (.NOT. (p.LE.0.5D0)) GO TO 30
      sign = -1.0D0
      z = p
      GO TO 40

   30 sign = 1.0D0
      z = 1.0D0 - dble(p)
   40 y = sqrt(-2.0D0*log(z))
      invnor = y + devlpl(xnum,5,y)/devlpl(xden,5,y)
      invnor = sign*invnor
      RETURN

      END
      LOGICAL FUNCTION LSAME ( CA, CB )
*     .. Scalar Arguments ..
      CHARACTER*1            CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME  tests if CA is the same letter as CB regardless of case.
*  CB is assumed to be an upper case letter. LSAME returns .TRUE. if
*  CA is either the same as CB or the equivalent lower case letter.
*
*  N.B. This version of the routine is only correct for ASCII code.
*       Installers must modify the routine for other character-codes.
*
*       For EBCDIC systems the constant IOFF must be changed to -64.
*       For CDC systems using 6-12 bit representations, the system-
*       specific code in comments must be activated.
*
*  Parameters
*  ==========
*
*  CA     - CHARACTER*1
*  CB     - CHARACTER*1
*           On entry, CA and CB specify characters to be compared.
*           Unchanged on exit.
*
*
*  Auxiliary routine for Level 2 Blas.
*
*  -- Written on 20-July-1986
*     Richard Hanson, Sandia National Labs.
*     Jeremy Du Croz, Nag Central Office.
*
*     .. Parameters ..
      INTEGER                IOFF
      PARAMETER            ( IOFF=32 )
*     .. Intrinsic Functions ..
      INTRINSIC              ICHAR
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA .EQ. CB
*
*     Now test for equivalence
*
      IF ( .NOT.LSAME ) THEN
         LSAME = ICHAR(CA) - IOFF .EQ. ICHAR(CB)
      END IF
*
      RETURN
*
*  The following comments contain code for CDC systems using 6-12 bit
*  representations.
*
*     .. Parameters ..
*     INTEGER                ICIRFX
*     PARAMETER            ( ICIRFX=62 )
*     .. Scalar Arguments ..
*     CHARACTER*1            CB
*     .. Array Arguments ..
*     CHARACTER*1            CA(*)
*     .. Local Scalars ..
*     INTEGER                IVAL
*     .. Intrinsic Functions ..
*     INTRINSIC              ICHAR, CHAR
*     .. Executable Statements ..
*
*     See if the first character in string CA equals string CB.
*
*     LSAME = CA(1) .EQ. CB .AND. CA(1) .NE. CHAR(ICIRFX)
*
*     IF (LSAME) RETURN
*
*     The characters are not identical. Now check them for equivalence.
*     Look for the 'escape' character, circumflex, followed by the
*     letter.
*
*     IVAL = ICHAR(CA(2))
*     IF (IVAL.GE.ICHAR('A') .AND. IVAL.LE.ICHAR('Z')) THEN
*        LSAME = CA(1) .EQ. CHAR(ICIRFX) .AND. CA(2) .EQ. CB
*     END IF
*
*     RETURN
*
*     End of LSAME.
*
      END
      real function rnor(jd)
c***begin prologue  rnor
c***date written   810915
c***revision date  830805
c***category no.  l6a14
c***keywords  random numbers, uniform random numbers
c***author    kahaner, david, scientific computing division, nbs
c             marsaglia, george, computer science dept., wash state univ
c
c***purpose  generates quasi normal random numbers, with mean zero and
c             unit standard deviation, and can be used with any computer
c             with integers at least as large as 32767.
c***description
c
c       rnor generates quasi normal random numbers with zero mean and
c       unit standard deviation.
c       it can be used with any computer with integers at least as
c       large as 32767.
c
c
c   use
c       first time....
c                   z = rnor(jd)
c                     here jd is any  n o n - z e r o  integer.
c                     this causes initialization of the program
c                     and the first random number to be returned as z.
c       subsequent times...
c                   z = rnor(0)
c                     causes the next random number to be returned as z.
c
c.....................................................................
c
c    note: users who wish to transport this program to other
c           computers should read the following ....
c
c   machine dependencies...
c      mdig = a lower bound on the number of binary digits available
c              for representing integers, including the sign bit.
c              this must be at least 16, but can be increased in
c              line with remark a below.
c
c   remarks...
c     a. this program can be used in two ways:
c        (1) to obtain repeatable results on different computers,
c            set 'mdig' to the smallest of its values on each, or,
c        (2) to allow the longest sequence of random numbers to be
c            generated without cycling (repeating) set 'mdig' to the
c            largest possible value.
c     b. the sequence of numbers generated depends on the initial
c          input 'jd' as well as the value of 'mdig'.
c          if mdig=16 one should find that
c            the first evaluation
c              z=rnor(87) gives  z=-.40079207...
c            the second evaluation
c              z=rnor(0) gives   z=-1.8728870...
c            the third evaluation
c              z=rnor(0) gives   z=1.8216004...
c            the fourth evaluation
c              z=rnor(0) gives   z=.69410355...
c            the thousandth evaluation
c              z=rnor(0) gives   z=.96782424...
c
c***references  marsaglia & tsang, "a fast, easily implemented
c                 method for sampling from decreasing or
c                 symmetric unimodal density functions", to be
c                 published in siam j sisc 1983.
c***routines called  i1mach,xerror
c***end prologue  rnor
      real v(65),w(65)
      integer m(17)
      save i1,j1,m,m1,m2,rmax
      data aa,b,c,rmax/12.37586,.4878992,12.67706,3.0518509e-5/
      data c1,c2,pc,xn/.9689279,1.301198,.1958303e-1,2.776994/
      data v/ .3409450, .4573146, .5397793, .6062427, .6631691
     +, .7136975, .7596125, .8020356, .8417227, .8792102, .9148948
     +, .9490791, .9820005, 1.0138492, 1.0447810, 1.0749254, 1.1043917
     +,1.1332738, 1.1616530, 1.1896010, 1.2171815, 1.2444516, 1.2714635
     +,1.2982650, 1.3249008, 1.3514125, 1.3778399, 1.4042211, 1.4305929
     +,1.4569915, 1.4834526, 1.5100121, 1.5367061, 1.5635712, 1.5906454
     +,1.6179680, 1.6455802, 1.6735255, 1.7018503, 1.7306045, 1.7598422
     +,1.7896223, 1.8200099, 1.8510770, 1.8829044, 1.9155830, 1.9492166
     +,1.9839239, 2.0198430, 2.0571356, 2.0959930, 2.1366450, 2.1793713
     +,2.2245175, 2.2725185, 2.3239338, 2.3795007, 2.4402218, 2.5075117
     +,2.5834658, 2.6713916, 2.7769943, 2.7769943, 2.7769943, 2.7769943/
      data w/   .10405134e-04, .13956560e-04, .16473259e-04,
     + .18501623e-04, .20238931e-04, .21780983e-04, .23182241e-04,
     + .24476931e-04, .25688121e-04, .26832186e-04, .27921226e-04,
     + .28964480e-04, .29969191e-04, .30941168e-04, .31885160e-04,
     + .32805121e-04, .33704388e-04, .34585827e-04, .35451919e-04,
     + .36304851e-04, .37146564e-04, .37978808e-04, .38803170e-04,
     + .39621114e-04, .40433997e-04, .41243096e-04, .42049621e-04,
     + .42854734e-04, .43659562e-04, .44465208e-04, .45272764e-04,
     + .46083321e-04, .46897980e-04, .47717864e-04, .48544128e-04,
     + .49377973e-04, .50220656e-04, .51073504e-04, .51937936e-04,
     + .52815471e-04, .53707761e-04, .54616606e-04, .55543990e-04,
     + .56492112e-04, .57463436e-04, .58460740e-04, .59487185e-04,
     + .60546402e-04, .61642600e-04, .62780711e-04, .63966581e-04,
     + .65207221e-04, .66511165e-04, .67888959e-04, .69353880e-04,
     + .70922996e-04, .72618816e-04, .74471933e-04, .76525519e-04,
     + .78843526e-04, .81526890e-04, .84749727e-04,
     + .84749727e-04, .84749727e-04, .84749727e-04/
      data m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),
     1     m(12),m(13),m(14),m(15),m(16),m(17)
     2                   / 30788,23052,2053,19346,10646,19427,23975,
     3                     19049,10949,19693,29746,26748,2796,23890,
     4                     29168,31924,16499 /
      data m1,m2,i1,j1 / 32767,256,5,17 /
c fast part...
c
c
c***first executable statement  rnor
      if(jd.ne.0)go to 27
   10 continue
      i=m(i1)-m(j1)
      if(i .lt. 0) i=i+m1
      m(j1)=i
      i1=i1-1
      if(i1 .eq. 0) i1=17
      j1=j1-1
      if(j1 .eq.. 0) j1=17
      j=mod(i,64)+1
      rnor=i*w(j+1)
      if( ( (i/m2)/2 )*2.eq.(i/m2))rnor=-rnor
      if(abs(rnor).le.v(j))return
c slow part; aa is a*f(0)
      x=(abs(rnor)-v(j))/(v(j+1)-v(j))
      y=uni(0)
      s=x+y
      if(s.gt.c2)go to 11
      if(s.le.c1)return
      if(y.gt.c-aa*exp(-.5*(b-b*x)**2))go to 11
      if(exp(-.5*v(j+1)**2)+y*pc/v(j+1).le.exp(-.5*rnor**2))return
c tail part; 3.855849 is .5*xn**2
   22 s=xn-alog(uni(0))/xn
      if(3.855849+alog(uni(0))-xn*s.gt.-.5*s**2)go to 22
      rnor=sign(s,rnor)
      return
   11 rnor=sign(b-b*x,rnor)
      return
c  fill
   27 continue
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      mdig=32
C     mdig=i1mach(8)+1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c          be sure that mdig at least 16...
C     if(mdig.lt.16)call xerror('rnor--mdig less than 16',23,1,2)
      m1 = 2**(mdig-2) + (2**(mdig-2)-1)
      m2 = 2**(mdig/2)
      jseed = min0(iabs(jd),m1)
      if( mod(jseed,2).eq.0 ) jseed=jseed-1
      k0 =mod(9069,m2)
      k1 = 9069/m2
      j0 = mod(jseed,m2)
      j1 = jseed/m2
      do 2 i=1,17
        jseed = j0*k0
        j1 = mod(jseed/m2+j0*k1+j1*k0,m2/2)
        j0 = mod(jseed,m2)
    2   m(i) = j0+m2*j1
      j1=17
      i1=5
      rmax = 1./float(m1)
c        seed uniform (0,1) generator.  (just a dummy call)
      rnor=uni(jd)
      do 28 i=1,65
   28  w(i)=rmax*v(i)
      go to 10
      end
      real function uni(jd)
c***begin prologue  uni
c***date written   810915
c***revision date  830805
c***category no.  l6a21
c***keywords  random numbers, uniform random numbers
c***author    blue, james, scientific computing division, nbs
c             kahaner, david, scientific computing division, nbs
c             marsaglia, george, computer science dept., wash state univ
c
c***purpose  this routine generates quasi uniform random numbers on [0,1
c             and can be used on any computer with which allows integers
c             at least as large as 32767.
c***description
c
c       this routine generates quasi uniform random numbers on the inter
c       [0,1).  it can be used with any computer which allows
c       integers at least as large as 32767.
c
c
c   use
c       first time....
c                   z = uni(jd)
c                     here jd is any  n o n - z e r o  integer.
c                     this causes initialization of the program
c                     and the first random number to be returned as z.
c       subsequent times...
c                   z = uni(0)
c                     causes the next random number to be returned as z.
c
c
c..................................................................
c   note: users who wish to transport this program from one computer
c         to another should read the following information.....
c
c   machine dependencies...
c      mdig = a lower bound on the number of binary digits available
c              for representing integers, including the sign bit.
c              this value must be at least 16, but may be increased
c              in line with remark a below.
c
c   remarks...
c     a. this program can be used in two ways:
c        (1) to obtain repeatable results on different computers,
c            set 'mdig' to the smallest of its values on each, or,
c        (2) to allow the longest sequence of random numbers to be
c            generated without cycling (repeating) set 'mdig' to the
c            largest possible value.
c     b. the sequence of numbers generated depends on the initial
c          input 'jd' as well as the value of 'mdig'.
c          if mdig=16 one should find that
c            the first evaluation
c              z=uni(305) gives z=.027832881...
c            the second evaluation
c              z=uni(0) gives   z=.56102176...
c            the third evaluation
c              z=uni(0) gives   z=.41456343...
c            the thousandth evaluation
c              z=uni(0) gives   z=.19797357...
c
c***references  marsaglia g., "comments on the perfect uniform random
c                 number generator", unpublished notes, wash s. u.
c***routines called  i1mach,xerror
c***end prologue  uni
      integer m(17)
c
      save i,j,m,m1,m2
c
      data m(1),m(2),m(3),m(4),m(5),m(6),m(7),m(8),m(9),m(10),m(11),
     1     m(12),m(13),m(14),m(15),m(16),m(17)
     2                   / 30788,23052,2053,19346,10646,19427,23975,
     3                     19049,10949,19693,29746,26748,2796,23890,
     4                     29168,31924,16499 /
      data m1,m2,i,j / 32767,256,5,17 /
c***first executable statement  uni
      if(jd .eq. 0) go to 3
c  fill
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      mdig=32
C     mdig=i1mach(8)+1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c          be sure that mdig at least 16...
C     if(mdig.lt.16)call xerror('uni--mdig less than 16',22,1,2)
      m1= 2**(mdig-2) + (2**(mdig-2)-1)
      m2 = 2**(mdig/2)
      jseed = min0(iabs(jd),m1)
      if( mod(jseed,2).eq.0 ) jseed=jseed-1
      k0 =mod(9069,m2)
      k1 = 9069/m2
      j0 = mod(jseed,m2)
      j1 = jseed/m2
      do 2 i=1,17
        jseed = j0*k0
        j1 = mod(jseed/m2+j0*k1+j1*k0,m2/2)
        j0 = mod(jseed,m2)
    2   m(i) = j0+m2*j1
      i=5
      j=17
c  begin main loop here
    3 k=m(i)-m(j)
      if(k .lt. 0) k=k+m1
      m(j)=k
      i=i-1
      if(i .eq. 0) i=17
      j=j-1
      if(j .eq. 0) j=17
      uni=float(k)/float(m1)
      return
      end
      SUBROUTINE XERBLA ( SRNAME, INFO )
*     ..    Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*6        SRNAME
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the Level 2 BLAS routines.
*
*  It is called by the Level 2 BLAS routines if an input parameter is
*  invalid.
*
*  Installers should consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Parameters
*  ==========
*
*  SRNAME - CHARACTER*6.
*           On entry, SRNAME specifies the name of the routine which
*           called XERBLA.
*
*  INFO   - INTEGER.
*           On entry, INFO specifies the position of the invalid
*           parameter in the parameter-list of the calling routine.
*
*
*  Auxiliary routine for Level 2 Blas.
*
*  Written on 20-July-1986.
*
*     .. Executable Statements ..
*
      WRITE (*,99999) SRNAME, INFO
*
      STOP
*
99999 FORMAT ( ' ** On entry to ', A6, ' parameter number ', I2,
     $         ' had an illegal value' )
*
*     End of XERBLA.
*
      END
*
