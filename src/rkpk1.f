      subroutine dbimdr (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y,
     * tol1, tol2, init, prec1, maxiter1, prec2, maxiter2, theta, nlaht,
     * score, varht, c, d, eta, wk, swk, qwk, ywk, u, w, info)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxiter1, 
     *maxiter2, info
      double precision s(lds,*), q(ldqr,ldqc,*), y(2,*), tol1, tol2, 
     *prec1, prec2, theta(*), nlaht, score(*), varht(2), c(*), d(*), 
     *wk(*), eta(*), swk(lds,*), qwk(ldqr,ldqc,*), ywk(*), u(*), w(*)
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
      varht(1) = 0.d0
      vmu=2
      j=1
23027 if(.not.(j.le.nobs))goto 23029
      varht(1) = varht(1) + u(j)**2 / w(j)
      j=j+1
      goto 23027
23029 continue
      varht(1) = varht(1) / dble (nobs)
23025 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dmudr (vmu, swk, lds, nobs, nnull, qwk, ldqr, ldqc, nq, ywk, 
     *tol2, init, prec1, maxiter1, theta, nlaht, score, varht, c, d, wk,
     * info)
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
cc    23005 goto 23004
      goto 23004
23006 continue
      return
      end
      subroutine dbisdr (vmu, s, lds, nobs, nnull, y, q, ldq, tol1, 
     *tol2, job, limnla, prec, maxiter, nlaht, score, varht, c, d, eta, 
     *qraux, jpvt, wk, swk, qwk, ywk, u, w, info)
c      character*2 vmu
      integer vmu
      integer lds, nobs, nnull, ldq, job, jpvt(*), info, maxiter
      double precision s(lds,*), y(2,*), q(ldq,*), tol1, tol2, limnla(2)
     *, nlaht, score(*), varht(2), c(*), d(*), qraux(*), wk(*), prec, et
     *a(*), swk(lds,*), qwk(ldq,*), ywk(*), u(*), w(*)
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
      varht(1) = 0.d0
      vmu=2
      j=1
23024 if(.not.(j.le.nobs))goto 23026
      varht(1) = varht(1) + u(j)**2 / w(j)
      j=j+1
      goto 23024
23026 continue
      varht(1) = varht(1) / dble (nobs)
23022 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dsidr (vmu, swk, lds, nobs, nnull, ywk, qwk, ldq, tol2, job, 
     *limnla, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
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
cc    23005 goto 23004
      goto 23004
23006 continue
      return
      end
      subroutine dbmdr (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, 
     *tol1, tol2, init, prec1, maxiter1, prec2, maxiter2, theta, nlaht, 
     *score, varht, c, d, eta, wk, swk, qwk, ywk, u, w, info)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxiter1, 
     *maxiter2, info
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol1, tol2, 
     *prec1, prec2, theta(*), nlaht, score(*), varht(2), c(*), d(*), 
     *wk(*), 
     *eta(*), swk(lds,*), qwk(ldqr,ldqc,*), ywk(*), u(*), w(*)
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
      varht(1) = 0.d0
      vmu=2
      j=1
23027 if(.not.(j.le.nobs))goto 23029
      varht(1) = varht(1) + u(j)**2 / w(j)
      j=j+1
      goto 23027
23029 continue
      varht(1) = varht(1) / dble (nobs)
23025 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dmudr (vmu, swk, lds, nobs, nnull, qwk, ldqr, ldqc, nq, ywk, 
     *tol2, init, prec1, maxiter1, theta, nlaht, score, varht, c, d, wk,
     * info)
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
cc    23005 goto 23004
      goto 23004
23006 continue
      return
      end
      subroutine dbsdr (vmu, s, lds, nobs, nnull, y, q, ldq, tol1, tol2,
     * job, limnla, prec, maxiter, nlaht, score, varht, c, d, eta, 
     *qraux, jpvt, wk, swk, qwk, ywk, u, w, info)
c      character*2 vmu
      integer vmu
      integer lds, nobs, nnull, ldq, job, jpvt(*), info, maxiter
      double precision s(lds,*), y(*), q(ldq,*), tol1, tol2, limnla(2), 
     *nlaht, score(*), varht(2), c(*), d(*), qraux(*), wk(*), prec, eta(
     **), swk(lds,*), qwk(ldq,*), ywk(*), u(*), w(*)
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
     *limnla, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
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
cc    23005 goto 23004
      goto 23004
23006 continue
      return
      end
      subroutine ddeev (vmu, nobs, q, ldqr, ldqc, n, nq, u, ldu, uaux, 
     *t, x, theta, nlaht, score, varht, hes, ldh, gra, hwk1, hwk2, gwk1,
     * gwk2, kwk, ldk, work1, work2, work3, info)
c      character*1 vmu
      integer vmu
      integer nobs, ldqr, ldqc, n, nq, ldu, ldh, ldk, info
      double precision q(ldqr,ldqc,*), u(ldu,*), uaux(*), t(2,*), x(*), 
     *theta(*), nlaht, score(*), varht(2), hes(ldh,*), gra(*), 
     *hwk1(nq,*), 
     *hwk2(nq,*), gwk1(*), gwk2(*), kwk(ldk,ldk,*), work1(*), work2(*), 
     *work3(*)
      double precision trc, det, dum, ddot, wk(1)
      integer i, j, m
      info = 0
      det = 0.d0
      dum = 0.d0
      trc = 0.d0
      call dset (nq, 0.d0, gra, 1)
      call dset (nq*nq, 0.d0, hes, 1)
c      if(.not.( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' ))
      if(.not.( vmu .ne. 0 .and. vmu .ne. 1 .and. vmu .ne. 2 ))
     *goto 23000
      info = -3
      return
23000 continue
      if(.not.( nobs .lt. n .or. ldqr .lt. n .or. ldqc .lt. n .or. nq 
     *.le. 0 .or. ldu .lt. n-1 .or. ldh .lt. nq .or. ldk .lt. n ))goto 2
     *3002
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
     *work1)
      call dqrsl (u, ldu, n-1, n-2, uaux, kwk(2,1,i), wk, kwk(2,1,i), 
     *wk, wk, wk, 01000, info)
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
     * 1)
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
     * 1)
      j=1
23065 if(.not.(j.le.i))goto 23067
      if(.not.( theta(j) .le. -25.d0 ))goto 23068
      goto 23066
23068 continue
      call dgemv ('n', n, n, 1.d0, kwk(1,1,j), n, work1, 1, 0.d0, work3,
     * 1)
      hwk1(i,j) = 2.d0 * ddot (n, work2, 1, work3, 1)
      call dgemv ('t', n, n, 1.d0, kwk(1,1,j), n, work1, 1, 0.d0, work3,
     * 1)
      hwk1(i,j) = hwk1(i,j) + 2.d0 * ddot (n, work2, 1, work3, 1)
23066 j=j+1
      goto 23065
23067 continue
      call dgemv ('t', n, n, 1.d0, kwk(1,1,i), n, work1, 1, 0.d0, work2,
     * 1)
      j=1
23070 if(.not.(j.le.i))goto 23072
      if(.not.( theta(j) .le. -25.d0 ))goto 23073
      goto 23071
23073 continue
      call dgemv ('n', n, n, 1.d0, kwk(1,1,j), n, work1, 1, 0.d0, work3,
     * 1)
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
     * 1)
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
      trc = dble (nobs) * 10.d0 ** (-nlaht) * varht(1) / score(1)
      i=1
23112 if(.not.(i.le.nq))goto 23114
      if(.not.( theta(i) .le. -25.d0 ))goto 23115
      goto 23113
23115 continue
      gra(i) = gwk1(i) / trc / trc - 2.d0 * score(1) * gwk2(i) / trc / 
     *dble(nobs)
23113 i=i+1
      goto 23112
23114 continue
      call dscal (nq, dble (nobs), gra, 1)
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
      call dscal (nq, 1.d0/dble (n), gra, 1)
23117 continue
c      if(.not.( vmu .eq. 'm' ))
      if(.not.( vmu .eq. 1 ))goto 23124
      det = 10.d0 ** (-nlaht) * varht(2) / score(1)
      i=1
23126 if(.not.(i.le.nq))goto 23128
      if(.not.( theta(i) .le. -25.d0 ))goto 23129
      goto 23127
23129 continue
      gra(i) = gwk1(i) / det - dble (nobs) / dble (n) * score(1) 
     * *gwk2(i)
23127 i=i+1
      goto 23126
23128 continue
      call dscal (nq, 1.d0 / dble (nobs), gra, 1)
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
     *** 3 - 2.d0 * gwk1(j) * gwk2(i) / trc ** 3 - 2.d0 * score(1) * 
     *hwk2(i,j) / trc / dble (nobs) + 6.d0 * score(1) * gwk2(i) * 
     *gwk2(j) / trc / trc / dble (nobs)
23139 j=j+1
      goto 23138
23140 continue
      call dscal (i, dble (nobs), hes(i,1), ldh)
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
     *j)
23151 j=j+1
      goto 23150
23152 continue
      call dscal (i, 1.d0/dble (n), hes(i,1), ldh)
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
      hes(i,j) = hwk1(i,j) / det - gwk1(i) * gwk2(j) / det / dble (n) 
     *- gwk1(j) * gwk2(i) / det / dble (n) - dble (nobs) / dble (
     *n) * score(1) * hwk2(i,j) + dble (nobs) / dble (n) ** 2 * 
     *score(1) *
     * gwk2(i) * gwk2(j)
23163 j=j+1
      goto 23162
23164 continue
      call dscal (i, 1.d0 / dble (nobs), hes(i,1), ldh)
23158 i=i+1
      goto 23157
23159 continue
23155 continue
      return
      end
      subroutine deval (vmu, q, ldq, n, z, nint, low, upp, nlaht, score,
     * varht, info, twk, work)
c      character*1 vmu
      integer vmu
      integer ldq, n, nint, info
      double precision q(ldq,*), z(*), low, upp, nlaht, score(*), varht(
     *2), twk(2,*), work(*)
      double precision tmp, minscr, mlo, varhtwk(2)
      integer j
      info = 0
      minscr = 0.d0
      varhtwk(1) = 0.d0
      varhtwk(2) = 0.d0
      if(.not.( upp .lt. low ))goto 23000
      mlo = low
      low = upp
      upp = mlo
23000 continue
c      if(.not.( (vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u')
      if(.not.( (vmu .ne. 0 .and. vmu .ne. 1 .and. vmu .ne. 2)  
     *.or. nint .lt. 1 ))goto 23002
      info = -3
      return
23002 continue
      if(.not.( 1 .gt. n .or. n .gt. ldq ))goto 23004
      info = -1
      return
23004 continue
      j=1
23006 if(.not.(j.le.nint+1))goto 23008
      tmp = low + dble (j-1) * ( upp - low ) / dble (nint)
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
      subroutine dgmdr (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, 
     *tol1, tol2, init, prec1, maxiter1, prec2, maxiter2, theta, nlaht, 
     *score, varht, c, d, eta, wk, swk, qwk, ywk, u, w, info)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxiter1, 
     *maxiter2, info
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol1, tol2, 
     *prec1, prec2, theta(*), nlaht, score(*), varht(2), c(*), d(*), 
     *wk(*), 
     *eta(*), swk(lds,*), qwk(ldqr,ldqc,*), ywk(*), u(*), w(*)
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
      varht(1) = 0.d0
      vmu=2
      j=1
23027 if(.not.(j.le.nobs))goto 23029
      varht(1) = varht(1) + u(j)**2 / w(j)
      j=j+1
      goto 23027
23029 continue
      varht(1) = varht(1) / dble (nobs)
23025 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dmudr (vmu, swk, lds, nobs, nnull, qwk, ldqr, ldqc, nq, ywk, 
     *tol2, init, prec1, maxiter1, theta, nlaht, score, varht, c, d, wk,
     * info)
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
cc    23005 goto 23004
      goto 23004
23006 continue
      return
      end
      subroutine dgold (vmu, q, ldq, n, z, low, upp, nlaht, score, 
     *varht, info, twk, work)
c      character*1 vmu
      integer vmu
      integer ldq, n, info
      double precision q(ldq,*), z(*), low, upp, nlaht, score(*), 
     * varht(2), twk(2,*), work(*)
      double precision ratio, mlo, mup, tmpl(1), tmpu(1)
      ratio = ( dsqrt (5.d0) - 1.d0 ) / 2.d0
      info = 0
      if(.not.( upp .lt. low ))goto 23000
      mlo = low
      low = upp
      upp = mlo
23000 continue
c      if(.not.( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' ))
      if(.not.( vmu .ne. 0 .and. vmu .ne. 1 .and. vmu .ne. 2 ))
     *goto 23002
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
      if(.not.( tmpl(1) .lt. tmpu(1) ))goto 23015
      upp = mup
      mup = mlo
      tmpu(1) = tmpl(1)
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
      tmpl(1) = tmpu(1)
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
cc    23011 goto 23010
      goto 23010
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
     * job, limnla, prec, maxiter, nlaht, score, varht, c, d, eta, 
     *qraux, jpvt, wk, swk, qwk, ywk, u, w, info)
c      character*2 vmu
      integer vmu
      integer lds, nobs, nnull, ldq, job, jpvt(*), info, maxiter
      double precision s(lds,*), y(*), q(ldq,*), tol1, tol2, limnla(2), 
     *nlaht, score(*), varht(2), c(*), d(*), qraux(*), wk(*), prec, 
     *eta(*),
     * swk(lds,*), qwk(ldq,*), ywk(*), u(*), w(*)
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
c      tmp = 1.d0
      tmp=dexp(-700.d0)
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
      varht(1) = 0.d0
      vmu=2
      j=1
23024 if(.not.(j.le.nobs))goto 23026
      varht(1) = varht(1) + u(j)**2 / w(j)
      j=j+1
      goto 23024
23026 continue
      varht(1) = varht(1) / dble (nobs)
23022 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dsidr (vmu, swk, lds, nobs, nnull, ywk, qwk, ldq, tol2, job, 
     *limnla, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
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
cc    23005 goto 23004
      goto 23004
23006 continue
      return
      end
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
      tmp = dsqrt (dble (p*p-1))
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
      delta = dasum (p, a, lda+1) / dble (p) * 1.d-7
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
     *tol, init, prec, maxite, theta, nlaht, score, varht, c, d, wk, 
     *info)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxite, 
     *jpvt(nnull), pvtwk(nnull), info
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol, prec, theta(
     **), nlaht, score(*), varht(2), c(*), d(*), wk(*)
c      character*1 vmu
      integer vmu
      integer n, n0
      integer iqraux, itraux, itwk, iqwk, iywk, ithewk, ihes, igra, 
     *ihwk1, ihwk2, igwk1, igwk2, ikwk, iwork1, iwork2
c    *ijpvt, ipvtwk
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
c      ijpvt = iwork2 + n
c      ipvtwk = ijpvt + n0
      call dmudr1 (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, tol, 
     *init, prec, maxite, theta, nlaht, score, varht, c, d, wk(iqraux), 
     *jpvt, wk(itwk), wk(itraux), wk(iqwk), wk(iywk), wk(ithewk), 
     *wk(ihes), wk(igra), wk(ihwk1), wk(ihwk2), wk(igwk1), wk(igwk2), 
     *pvtwk, wk(ikwk), wk(iwork1), wk(iwork2), info)
c      call dmudr1 (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y, tol, 
c     *init, prec, maxite, theta, nlaht, score, varht, c, d, wk(iqraux), 
c     *wk(ijpvt), wk(itwk), wk(itraux), wk(iqwk), wk(iywk), wk(ithewk), 
c     *wk(ihes), wk(igra), wk(ihwk1), wk(ihwk2), wk(igwk1), wk(igwk2), 
c     *wk(ipvtwk), wk(ikwk), wk(iwork1), wk(iwork2), info)
      return
      end
      subroutine dmudrnew (vmu, s, lds, nobs, nnull, q,q1, q2, ldqr, 
     *ldqc,nq, y, tol, init, prec, maxite, theta, nlaht, score, varht, 
     *c, d, wk, ok, info)
      integer lds, nobs, nnull, ldqr, ldqc, init, maxite, info, nq, ok
      double precision s(lds,*),q(ldqr,ldqc,*), y(*), tol, prec, theta(*
     *), nlaht, score(*), varht(2), c(*), d(*), wk(*), q1(ldqr,*), 
     *q2(ldqr,*)
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
     *init, prec, maxite, theta, nlaht, score, varht, c, d, wk, info)
      if(.not.(d(1).lt. 0.001))goto 23006
      ok=1
23006 continue
      return
      end
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
     *tol1, tol2, init, prec1, maxiter1, prec2, maxiter2, theta, nlaht, 
     *score, varht, c, d, eta, wk, swk, qwk, ywk, u, w, info)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxiter1, 
     *maxiter2, info
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol1, tol2, 
     *prec1, prec2, theta(*), nlaht, score(*), varht(2), c(*), d(*), 
     *wk(*), 
     *eta(*), swk(lds,*), qwk(ldqr,ldqc,*), ywk(*), u(*), w(*)
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
      varht(1) = 0.d0
      vmu=2
      j=1
23027 if(.not.(j.le.nobs))goto 23029
      varht(1) = varht(1) + u(j)**2 / w(j)
      j=j+1
      goto 23027
23029 continue
      varht(1) = varht(1) / dble (nobs)
23025 continue
      call dcopy (nobs, ywk, 1, u, 1)
      call dmudr (vmu, swk, lds, nobs, nnull, qwk, ldqr, ldqc, nq, ywk, 
     *tol2, init, prec1, maxiter1, theta, nlaht, score, varht, c, d, wk,
     * info)
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
cc    23005 goto 23004
      goto 23004
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
     * job, limnla, prec, maxiter, nlaht, score, varht, c, d, eta, 
     *qraux, jpvt, wk, swk, qwk, ywk, u, w, info)
c      character*2 vmu
      integer vmu
      integer lds, nobs, nnull, ldq, job, jpvt(*), info, maxiter
      double precision s(lds,*), y(*), q(ldq,*), tol1, tol2, limnla(2), 
     *nlaht, score(*), varht(2), c(*), d(*), qraux(*), wk(*), prec, eta(
     **), swk(lds,*), qwk(ldq,*), ywk(*), u(*), w(*)
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
     *limnla, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
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
cc    23005 goto 23004
      goto 23004
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
c     dqrsl applies the output of dqrdc to compute coordinate
c     transformations, projections, and least squares solutions.
c     for k .le. min(n,p), let xk be the matrix
c
c            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))
c
c     formed from columnns jpvt(1), ... ,jpvt(k) of the original
c     n x p matrix x that was input to dqrdc (if no pivoting was
c     done, xk consists of the first k columns of x in their
c     original order).  dqrdc produces a factored orthogonal matrix q
c     and an upper triangular matrix r such that
c
c              xk = q * (r)
c                       (0)
c
c     this information is contained in coded form in the arrays
c     x and qraux.
c
c     on entry
c
c        x      double precision(ldx,p).
c               x contains the output of dqrdc.
c
c        ldx    integer.
c               ldx is the leading dimension of the array x.
c
c        n      integer.
c               n is the number of rows of the matrix xk.  it must
c               have the same value as n in dqrdc.
c
c        k      integer.
c               k is the number of columns of the matrix xk.  k
c               must nnot be greater than min(n,p), where p is the
c               same as in the calling sequence to dqrdc.
c
c        qraux  double precision(p).
c               qraux contains the auxiliary output from dqrdc.
c
c        y      double precision(n)
c               y contains an n-vector that is to be manipulated
c               by dqrsl.
c
c        job    integer.
c               job specifies what is to be computed.  job has
c               the decimal expansion abcde, with the following
c               meaning.
c
c                    if a.ne.0, compute qy.
c                    if b,c,d, or e .ne. 0, compute qty.
c                    if c.ne.0, compute b.
c                    if d.ne.0, compute rsd.
c                    if e.ne.0, compute xb.
c
c               note that a request to compute b, rsd, or xb
c               automatically triggers the computation of qty, for
c               which an array must be provided in the calling
c               sequence.
c
c     on return
c
c        qy     double precision(n).
c               qy conntains q*y, if its computation has been
c               requested.
c
c        qty    double precision(n).
c               qty contains trans(q)*y, if its computation has
c               been requested.  here trans(q) is the
c               transpose of the matrix q.
c
c        b      double precision(k)
c               b contains the solution of the least squares problem
c
c                    minimize norm2(y - xk*b),
c
c               if its computation has been requested.  (note that
c               if pivoting was requested in dqrdc, the j-th
c               component of b will be associated with column jpvt(j)
c               of the original matrix x that was input into dqrdc.)
c
c        rsd    double precision(n).
c               rsd contains the least squares residual y - xk*b,
c               if its computation has been requested.  rsd is
c               also the orthogonal projection of y onto the
c               orthogonal complement of the column space of xk.
c
c        xb     double precision(n).
c               xb contains the least squares approximation xk*b,
c               if its computation has been requested.  xb is also
c               the orthogonal projection of y onto the column space
c               of x.
c
c        info   integer.
c               info is zero unless the computation of b has
c               been requested and r is exactly singular.  in
c               this case, info is the index of the first zero
c               diagonal element of r and b is left unaltered.
c
c     the parameters qy, qty, b, rsd, and xb are not referenced
c     if their computation is not requested and in this case
c     can be replaced by dummy variables in the calling program.
c     to save storage, the user may in some cases use the same
c     array for different parameters in the calling sequence.  a
c     frequently occuring example is when one wishes to compute
c     any of b, rsd, or xb and does not need y or qty.  in this
c     case one may identify y, qty, and one of b, rsd, or xb, while
c     providing separate arrays for anything else that is to be
c     computed.  thus the calling sequence
c
c          call dqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)
c
c     will result in the computation of b and rsd, with rsd
c     overwriting y.  more generally, each item in the following
c     list contains groups of permissible identifications for
c     a single callinng sequence.
c
c          1. (y,qty,b) (rsd) (xb) (qy)
c
c          2. (y,qty,rsd) (b) (xb) (qy)
c
c          3. (y,qty,xb) (b) (rsd) (qy)
c
c          4. (y,qy) (qty,b) (rsd) (xb)
c
c          5. (y,qy) (qty,rsd) (b) (xb)
c
c          6. (y,qy) (qty,xb) (b) (rsd)
c
c     in any group the value returned in the array allocated to
c     the group corresponds to the last member of the group.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dqrsl uses the following functions and subprograms.
c
c     BLAS      daxpy,dcopy,ddot
c     Fortran   dabs,min0,mod
c
      subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
      integer ldx,n,k,job,info
      double precision x(ldx,*),qraux(*),y(*),qy(*),qty(*),b(*),rsd(*),
     *                 xb(*)
c
c     internal variables
c
      integer i,j,jj,ju,kp1
      double precision ddot,t,temp
      logical cb,cqy,cqty,cr,cxb
c
c
c     set info flag.
c
      info = 0
c
c     determine what is to be computed.
c
      cqy = job/10000 .ne. 0
      cqty = mod(job,10000) .ne. 0
      cb = mod(job,1000)/100 .ne. 0
      cr = mod(job,100)/10 .ne. 0
      cxb = mod(job,10) .ne. 0
      ju = min0(k,n-1)
c
c     special action when n=1.
c
      if (ju .ne. 0) go to 40
         if (cqy) qy(1) = y(1)
         if (cqty) qty(1) = y(1)
         if (cxb) xb(1) = y(1)
         if (.not.cb) go to 30
            if (x(1,1) .ne. 0.0d0) go to 10
               info = 1
            go to 20
 10                continue
               b(1) = y(1)/x(1,1)
 20                   continue
 30                       continue
         if (cr) rsd(1) = 0.0d0
      go to 250
 40    continue
c
c        set up to compute qy or qty.
c
         if (cqy) call dcopy(n,y,1,qy,1)
         if (cqty) call dcopy(n,y,1,qty,1)
         if (.not.cqy) go to 70
c
c           compute qy.
c
            do 60 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 50
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
                  call daxpy(n-j+1,t,x(j,j),1,qy(j),1)
                  x(j,j) = temp
 50                         continue
 60                                continue
 70                                    continue
         if (.not.cqty) go to 100
c
c           compute trans(q)*y.
c
            do 90 j = 1, ju
               if (qraux(j) .eq. 0.0d0) go to 80
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
                  call daxpy(n-j+1,t,x(j,j),1,qty(j),1)
                  x(j,j) = temp
 80                         continue
 90                                continue
 100                                   continue
c
c        set up to compute b, rsd, or xb.
c
         if (cb) call dcopy(k,qty,1,b,1)
         kp1 = k + 1
         if (cxb) call dcopy(k,qty,1,xb,1)
         if (cr .and. k .lt. n) call dcopy(n-k,qty(kp1),1,rsd(kp1),1)
         if (.not.cxb .or. kp1 .gt. n) go to 120
            do 110 i = kp1, n
               xb(i) = 0.0d0
 110                  continue
 120                      continue
         if (.not.cr) go to 140
            do 130 i = 1, k
               rsd(i) = 0.0d0
 130                  continue
 140                      continue
         if (.not.cb) go to 190
c
c           compute b.
c
            do 170 jj = 1, k
               j = k - jj + 1
               if (x(j,j) .ne. 0.0d0) go to 150
                  info = j
c           ......exit
                  go to 180
 150                        continue
               b(j) = b(j)/x(j,j)
               if (j .eq. 1) go to 160
                  t = -b(j)
                  call daxpy(j-1,t,x(1,j),1,b,1)
 160                        continue
 170                               continue
 180                                      continue
 190                                          continue
         if (.not.cr .and. .not.cxb) go to 240
c
c           compute rsd or xb as required.
c
            do 230 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 220
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  if (.not.cr) go to 200
                     t = -ddot(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
                     call daxpy(n-j+1,t,x(j,j),1,rsd(j),1)
 200                              continue
                  if (.not.cxb) go to 210
                     t = -ddot(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
                     call daxpy(n-j+1,t,x(j,j),1,xb(j),1)
 210                              continue
                  x(j,j) = temp
 220                        continue
 230                               continue
 240                                   continue
 250                                    continue
      return
      end
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
      double precision dum(1), ddot
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
      subroutine dstup (s, lds, nobs, nnull, qraux, jpvt, y, q, ldqr, 
     *ldqc, nq, info, work, dum)
      integer lds, nobs, nnull, jpvt(*), ldqr, ldqc, nq, info
      double precision s(lds,*), y(*), qraux(*), q(ldqr,ldqc,*), work(*)
      double precision dum(1)
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
      subroutine dtrev (vmu, t, ldt, n, z, score, varht, info, work)      
c      character*1 vmu
      integer vmu
      integer n, info
      double precision t(ldt,*), z(*), score(*), varht(2), work(*)
      double precision nume, deno, tmp, alph, la, dasum, ddot
      integer j
      info = 0
c      if(.not.( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' ))
      if(.not.( vmu .ne. 0 .and. vmu .ne. 1 .and. vmu .ne. 2 ))
     *goto 23000
      info = -3
      return
23000 continue
      la = t(1,1)
      alph = dble (n) / dasum (n, t(2,1), ldt)
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
      nume = ddot (n, work, 1, work, 1) / dble (n)
      deno = deno / dble (n)
      varht(1) = alph * la * nume / deno
      score(1) = nume / deno / deno
      deno = dlog (t(2,n))
      j=n-1
23009 if(.not.(j.gt.0))goto 23011
      deno = deno + dlog (t(2,j))
      j=j-1
      goto 23009
23011 continue
      nume = ddot (n, z, 1, work, 1) / dble (n)
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
      nume = ddot (n, z, 1, work, 1) / dble (n)
      varht(2) = alph * la * nume
      score(1) = nume * dexp (2.d0 * deno / dble (n))
      tmp = 1.d0 / t(2,n) / t(2,n)
      deno = tmp
      j=n-1
23017 if(.not.(j.gt.0))goto 23019
      tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
      deno = deno + tmp
      j=j-1
      goto 23017
23019 continue
      nume = ddot (n, work, 1, work, 1) / dble (n)
      deno = deno / dble (n)
      varht(1) = alph * la * nume / deno
23012 continue
c      if(.not.( vmu .eq. 'u' ))
      if(.not.( vmu .eq. 2 ))goto 23020
      nume = ddot (n, work, 1, work, 1) / dble (n)
      tmp = 1.d0 / t(2,n) / t(2,n)
      deno = tmp
      j=n-1
23022 if(.not.(j.gt.0))goto 23024
      tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
      deno = deno + tmp
      j=j-1
      goto 23022
23024 continue
      deno = deno / dble (n)
      score(1) = alph * alph * la * la * nume - 2.d0 * varht(1) * 
     *alph * la * deno
      varht(2) = alph * la *deno
23020 continue
      return
      end

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
      dn = dble (n)
      toluni = toltot * 6.d0 / dn / ( dn - 1.d0 ) / ( 2.d0 * dn - 1.d0 )
      tolcum = 0.d0
      j=1 
23009 if(.not.(j.lt.n-1))goto 23011
      nrmtot = nrmtot - x(j,j) * x(j,j)
      nrmxj = ddot (n-j, x(j+1,j), 1, x(j+1,j), 1)
      dn = dble (n-j)
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
