      subroutine dcoef (s, lds, nobs, nnull, qraux, jpvt, z, q, ldq, 
&     nlaht, c, d, info, twk)
      integer lds, nobs, nnull, jpvt(*), ldq, info
      double precision s(lds,*), qraux(*), z(*), q(ldq,*), nlaht, c(*), 
&     d(*), twk(2,*)
      double precision dum, ddot
      integer n, n0
      info = 0
      if(.not.( nnull .lt. 0 .or. nnull .ge. nobs .or. nobs .gt. lds 
&     .or. nobs .gt. ldq ))goto 23000
      info = -1
      return
23000 continue
      n0 = nnull
      n = nobs - nnull
      call dset (n, 10.d0 ** nlaht, twk(2,1), 2)
      call daxpy (n, 1.d0, q(n0+1,n0+1), ldq+1, twk(2,1), 2)
      call dcopy (n-1, q(n0+1,n0+2), ldq+1, twk(1,2), 2)
      call dpbfa (twk, 2, n, 1, info)
      if(.not.( info .ne. 0 ))goto 23002
      info = -2
      return
23002 continue
      call dpbsl (twk, 2, n, 1, z(n0+1))
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, twk, 1)
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, twk, z(n0+2), z(n0+2), 
&     dum, dum, dum, dum, 10000, info)
      if(.not.( nnull .eq. 0 ))goto 23004
      call dcopy (n, z(n0+1), 1, c(n0+1), 1)
      return
23004 continue
      call dset (n0, 0.d0, c, 1)
      call dcopy (n, z(n0+1), 1, c(n0+1), 1)
      call dqrsl (s, lds, nobs, nnull, qraux, c, c, dum, dum, dum, dum, 
&     10000, info)
      j=1
23006 if(.not.(j.le.n0))goto 23008
      d(j) = z(j) - ddot (n, z(n0+1), 1, q(n0+1,j), 1)
      j=j+1
      goto 23006
23008 continue
      call dtrsl (s, lds, n0, d, 01, info)
      call dprmut (d, n0, jpvt, 1)
      return
      end
      subroutine dsidr (vmu, s, lds, nobs, nnull, y, q, ldq, tol, job, 
&     limnla, nlaht, score, varht, c, d, qraux, jpvt, wk, info)
c      character*1 vmu
      integer vmu
      integer lds, nobs, nnull, ldq, job, jpvt(*), info
      double precision s(lds,*), y(*), q(ldq,*), tol, limnla(2), nlaht, 
&     score(*), varht(2), c(*), d(*), qraux(*), wk(*), tmp
      info = 0
      if(.not.( nnull .lt. 0 .or. nnull .ge. nobs .or. nobs .gt. lds 
&     .or. nobs .gt. ldq ))goto 23000
      info = -1
      return
23000 continue
c      if(.not.( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' ))
      if(.not.( vmu .ne. 0 .and. vmu .ne. 1 .and. vmu .ne. 2 ))
&     goto 23002
      info = -3
      return
23002 continue
      if(.not.( nnull .eq. 0 ))goto 23004
      call dcore (vmu, q, ldq, nobs, nnull, tol, y, job, limnla, nlaht, 
&     score, varht, info, wk, wk(2*nobs+1))
      if(.not.( info .ne. 0 ))goto 23006
      return
23006 continue
      call dcoef (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nlaht, c,
&      d, info, wk)
      return
23004 continue
      call dstup (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nobs, 1, 
&     info, wk)
      if(.not.( info .ne. 0 ))goto 23008
      return
23008 continue
      call dcore (vmu, q, ldq, nobs, nnull, tol, y, job, limnla, nlaht, 
&     score, varht, info, wk, wk(2*nobs+1))
      if(.not.( info .ne. 0 ))goto 23010
      return
23010 continue
      call dcoef (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nlaht, c,
&      d, info, wk)
      return
      end
      subroutine dcore (vmu, q, ldq, nobs, nnull, tol, z, job, limnla, 
&     nlaht, score, varht, info, twk, work)
c      character*1 vmu
      integer vmu
      integer ldq, nobs, nnull, job, info
      double precision q(ldq,*), tol, z(*), limnla(2), nlaht, score(*), 
&     varht(2), twk(2,*), work(*)
      double precision dum, low, upp, dasum, mchpr
      integer n0, n, j
      info = 0
c      if(.not.( vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u' ))
      if(.not.( vmu .ne. 0 .and. vmu .ne. 1 .and. vmu .ne. 2 ))
&     goto 23000
      info = -3
      return
23000 continue
      if(.not.( nnull .lt. 0 .or. nobs .le. nnull .or. nobs .gt. ldq ))
&     goto 23002
      info = -1
      return
23002 continue
      n0 = nnull
      n = nobs - nnull
      call dsytr (q(n0+1,n0+1), ldq, n, tol, info, work)
      if(.not.( info .ne. 0 ))goto 23004
      return
23004 continue
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, work, 1)
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, work, z(n0+2), dum, z(n0+
&     2), dum, dum, dum, 01000, info)
      if(.not.( job .eq. 0 ))goto 23006
      mchpr = 1.d0
23008 if(.not.( 1.d0 + mchpr .gt. 1.d0 ))goto 23009
      mchpr = mchpr / 2.d0
      goto 23008
23009 continue
      mchpr = mchpr * 2.d0
      limnla(2) = dmax1 (dasum (n, q(n0+1,n0+1), ldq+1) * 1.d2, mchpr)
      limnla(1) = limnla(2) * mchpr
      limnla(2) = dlog10 (limnla(2))
      limnla(1) = dlog10 (limnla(1))
23006 continue
      low = limnla(1)
      upp = limnla(2)
      if(.not.( job .le. 0 ))goto 23010
      call dgold (vmu, q(n0+1,n0+1), ldq, n, z(n0+1), low, upp, nlaht, 
&     score(1), varht, info, twk, work)
c      if(.not.( vmu .eq. 'v' ))
      if(.not.( vmu .eq. 0 ))goto 23012
      score(1) = score(1) * dfloat (nobs) / dfloat (n)
23012 continue
c      if(.not.( vmu .eq. 'm' ))
      if(.not.( vmu .eq. 1 ))goto 23014  
      score(1) = score(1) * dfloat (n) / dfloat (nobs)
23014 continue
c      if(.not.( vmu .eq. 'u' ))
      if(.not.( vmu .eq. 2 ))goto 23016
      score(1) = score(1) * dfloat (n) / dfloat (nobs) + 2.d0 * varht(1)
23016 continue
      goto 23011
23010 continue
      call deval (vmu, q(n0+1,n0+1), ldq, n, z(n0+1), job, low, upp, 
&     nlaht, score, varht, info, twk, work)
      dum = dfloat (nobs) / dfloat (n)
      j=1
23018 if(.not.(j.le.job+1))goto 23020
c      if(.not.( vmu .eq. 'v' ))
      if(.not.( vmu .eq. 0 ))goto 23021
      score(j) = score(j) * dum
23021 continue
c      if(.not.( vmu .eq. 'm' ))
      if(.not.( vmu .eq. 1 ))goto 23023
      score(j) = score(j) / dum
23023 continue
c      if(.not.( vmu .eq. 'u' ))
      if(.not.( vmu .eq. 2 ))goto 23025
      score(j) = score(j) / dum + 2.d0 * varht(1)
23025 continue
      j=j+1
      goto 23018
23020 continue
23011 continue
      return
      end
      subroutine dmudr1 (vmu, s, lds, nobs, nnull, q, ldqr, ldqc, nq, y,
&      tol, init, prec, maxite, theta, nlaht, score, varht, c, d, qraux,
&      jpvt, twk, traux, qwk, ywk, thewk, hes, gra, hwk1, hwk2, gwk1, 
&     gwk2, pvtwk, kwk, work1, work2, info)
      integer lds, nobs, nnull, ldqr, ldqc, nq, init, maxite, jpvt(*), 
&     pvtwk(*), info
      double precision s(lds,*), q(ldqr,ldqc,*), y(*), tol, prec, theta(
&     *), nlaht, score, varht(2), c(*), d(*), qraux(*), traux(*), twk(2,
&     *), qwk(nobs,*), ywk(*), thewk(*), hes(nq,*), gra(*), hwk1(nq,*), 
&     hwk2(nq,*), gwk1(*), gwk2(*), kwk(nobs-nnull,nobs-nnull,*), work1(
&     *), work2(*)
c      character*1 vmu
      integer vmu
      double precision alph, scrold, scrwk, nlawk, limnla(2), tmp, 
&     dasum, ddot
      integer n, n0, i, j, iwk, maxitwk, idamax, job
      info = 0
      n0 = nnull
      n = nobs - nnull
      maxitwk = maxite
c      if(.not.( (vmu .ne. 'v' .and. vmu .ne. 'm' .and. vmu .ne. 'u') 
      if(.not.( (vmu .ne.  0  .and. vmu .ne.  1  .and. vmu .ne.  2 ) 
&     .or. (init .ne. 0 .and. init .ne. 1) .or. (maxitwk .le.0) .or. (
&     prec .le. 0.d0) ))goto 23000
      info = -3
      return
23000 continue
      if(.not.( lds .lt. nobs .or. nobs .le. n0 .or. n0 .lt. 0 .or. 
&     ldqr .lt. nobs .or. ldqc .lt. nobs .or. nq .le. 0 ))goto 23002
      info = -1
      return
23002 continue
      if(.not.( n0 .gt. 0 ))goto 23004
      call dstup (s, lds, nobs, n0, qraux, jpvt, y, q, ldqr, ldqc, nq, 
&     info, work1)
23004 continue
      if(.not.( info .ne. 0 ))goto 23006
      return
23006 continue
      if(.not.( init .eq. 1 ))goto 23008
      call dcopy (nq, theta, 1, thewk, 1)
      goto 23009
23008 continue
      i=1
23010 if(.not.(i.le.nq))goto 23012
      thewk(i) = dasum (n, q(n0+1,n0+1,i), ldqr+1)
      if(.not.( thewk(i) .gt. 0.d0 ))goto 23013
      thewk(i) = 1.d0 / thewk(i)
23013 continue
      i=i+1
      goto 23010
23012 continue
      j=1
23015 if(.not.(j.le.nobs))goto 23017
      call dset (nobs-j+1, 0.d0, qwk(j,j), 1)
      j=j+1
      goto 23015
23017 continue
      i=1
23018 if(.not.(i.le.nq))goto 23020
      j=1
23021 if(.not.(j.le.nobs))goto 23023
      call daxpy (nobs-j+1, thewk(i), q(j,j,i), 1, qwk(j,j), 1)
      j=j+1
      goto 23021
23023 continue
      i=i+1
      goto 23018
23020 continue
      call dcopy (nobs, y, 1, ywk, 1)
      call dcore (vmu, qwk, nobs, nobs, n0, tol, ywk, 0, limnla, nlawk, 
&     scrwk, varht, info, twk, work1)
      if(.not.(info .ne. 0 ))goto 23024
      return
23024 continue
      call dcoef (s, lds, nobs, n0, qraux, jpvt, ywk, qwk, nobs, nlawk, 
&     c, d, info, twk)
      if(.not.( n0 .gt. 0 ))goto 23026
      call dqrsl (s, lds, nobs, n0, qraux, c, tmp, c, tmp, tmp, tmp, 
&     01000, info)
23026 continue
      i=1
23028 if(.not.(i.le.nq))goto 23030
      call dsymv('l', n, thewk(i), q(n0+1,n0+1,i), ldqr, c(n0+1), 1, 0.
&     d0, work1, 1)
      thewk(i) = ddot (n, c(n0+1), 1, work1, 1) * thewk(i)
      if(.not.( thewk(i) .gt. 0.d0 ))goto 23031
      thewk(i) = dlog10 (thewk(i))
      goto 23032
23031 continue
      thewk(i) = -25.d0
23032 continue
      i=i+1
      goto 23028
23030 continue
23009 continue
      scrold = 1.d10
      job = 0
23033 continue
      if(.not.( nq .eq. 1 ))goto 23036
      theta(1) = 0.d0
      goto 23035
23036 continue
      j=1
23038 if(.not.(j.le.nobs))goto 23040
      call dset (nobs-j+1, 0.d0, qwk(j,j), 1)
      j=j+1
      goto 23038
23040 continue
      i=1
23041 if(.not.(i.le.nq))goto 23043
      if(.not.( thewk(i) .le. -25.d0 ))goto 23044
      goto 23042
23044 continue
      j=1
23046 if(.not.(j.le.nobs))goto 23048
      call daxpy (nobs-j+1, 10.d0 ** thewk(i), q(j,j,i), 1, qwk(j,j), 1)
      j=j+1
      goto 23046
23048 continue
23042 i=i+1
      goto 23041
23043 continue
      call dcopy (nobs, y, 1, ywk, 1)
      call dcore (vmu, qwk, nobs, nobs, n0, tol, ywk, job, limnla, 
&     nlawk, scrwk, varht, info, twk, work1)
      if(.not.(info .ne. 0 ))goto 23049
      return
23049 continue
      if(.not.( scrold .lt. scrwk ))goto 23051
      tmp = dabs (gwk1(idamax (nq, gwk1, 1)))
      if(.not.( alph * tmp .gt. - prec ))goto 23053
      info = -5
      return
23053 continue
      alph = alph / 2.d0
      i=1
23055 if(.not.(i.le.nq))goto 23057
      thewk(i) = theta(i) + alph * gwk1(i)
      i=i+1
      goto 23055
23057 continue
      goto 23034
23051 continue
      maxitwk = maxitwk - 1
      call dcopy (n-2, qwk(n0+2,n0+1), nobs+1, traux, 1)
      call dcopy (n, qwk(n0+1,n0+1), nobs+1, twk(2,1), 2)
      call dcopy (n-1, qwk(n0+1,n0+2), nobs+1, twk(1,2), 2)
      call ddeev (vmu, nobs, q(n0+1,n0+1,1), ldqr, ldqc, n, nq, qwk(n0+
&     2,n0+1), nobs, traux, twk, ywk(n0+1), thewk, nlawk, scrwk, varht, 
&     hes, nq, gra, hwk1, hwk2, gwk1, gwk2, kwk, n, work1, work2, c, 
&     info)
      iwk = 0
      i=1
23058 if(.not.(i.le.nq))goto 23060
      if(.not.( thewk(i) .le. -25.d0 ))goto 23061
      goto 23059
23061 continue
      iwk = iwk + 1
      call dcopy (nq, hes(1,i), 1, hes(1,iwk), 1)
23059 i=i+1
      goto 23058
23060 continue
      iwk = 0
      i=1
23063 if(.not.(i.le.nq))goto 23065
      if(.not.( thewk(i) .le. -25.d0 ))goto 23066
      goto 23064
23066 continue
      iwk = iwk + 1
      call dcopy (nq, hes(i,1), nq, hes(iwk,1), nq)
      gwk1(iwk) = gra(i)
      work2(iwk) = gra(i)
23064 i=i+1
      goto 23063
23065 continue
      i=1
23068 if(.not.(i.lt.iwk))goto 23070
      call dcopy (iwk-i, hes(i+1,i), 1, hes(i,i+1), nq)
      i=i+1
      goto 23068
23070 continue
      call dmcdc (hes, nq, iwk, gwk2, pvtwk, info)
      call dprmut (gwk1, iwk, pvtwk, 0)
      call dposl (hes, nq, iwk, gwk1)
      call dprmut (gwk1, iwk, pvtwk, 1)
      alph = -1.d0
      j = iwk
      i=nq
23071 if(.not.(i.ge.1))goto 23073
      if(.not.( thewk(i) .le. -25.0 ))goto 23074
      gwk1(i) = 0.d0
      goto 23075
23074 continue
      gwk1(i) = gwk1(iwk)
      iwk = iwk - 1
23075 continue
      i=i-1
      goto 23071
23073 continue
      call dscal (nq, 1.d0/dlog(1.d1), gwk1, 1)
      tmp = dabs (gwk1(idamax (nq, gwk1, 1)))
      if(.not.( tmp .gt. 1.d0 ))goto 23076
      call dscal (nq, 1.d0/tmp, gwk1, 1)
23076 continue
      i=1
23078 if(.not.(i.le.nq))goto 23080
      if(.not.( thewk(i) .le. -25.d0 ))goto 23081
      goto 23079
23081 continue
      thewk(i) = thewk(i) - nlawk
23079 i=i+1
      goto 23078
23080 continue
      call dcopy (nq, thewk, 1, theta, 1)
      tmp = gra(idamax (nq, gra, 1)) ** 2
      if(.not.( tmp .lt. prec ** 2 .or. scrold - scrwk .lt. prec * (
&     scrwk + 1.d0) .and. tmp .lt. prec * (scrwk + 1.d0) ** 2 ))goto 230
&     83
      goto 23035
23083 continue
      if(.not.( maxitwk .lt. 1 ))goto 23085
      info = -4
      return
23085 continue
      scrold = scrwk
      i=1
23087 if(.not.(i.le.nq))goto 23089
      thewk(i) = thewk(i) + alph * gwk1(i)
      i=i+1
      goto 23087
23089 continue
      job = -1
      limnla(1) = -1.d0
      limnla(2) = 1.d0
23034 goto 23033
23035 continue
      j=1
23090 if(.not.(j.le.nobs))goto 23092
      call dset (nobs-j+1, 0.d0, qwk(j,j), 1)
      j=j+1
      goto 23090
23092 continue
      i=1
23093 if(.not.(i.le.nq))goto 23095
      if(.not.( theta(i) .le. -25.d0 ))goto 23096
      goto 23094
23096 continue
      j=1
23098 if(.not.(j.le.nobs))goto 23100
      call daxpy (nobs-j+1, 10.d0 ** theta(i), q(j,j,i), 1, qwk(j,j), 1)
      j=j+1
      goto 23098
23100 continue
23094 i=i+1
      goto 23093
23095 continue
      call dcopy (nobs, y, 1, ywk, 1)
      call dcore (vmu, qwk, nobs, nobs, n0, tol, ywk, job, limnla, 
&     nlaht, score, varht, info, twk, work1)
      if(.not.(info .ne. 0 ))goto 23101
      return
23101 continue
      call dcoef (s, lds, nobs, n0, qraux, jpvt, ywk, qwk, nobs, nlaht, 
&     c, d, info, twk)
      return
      end
      subroutine dcrdr (s, lds, nobs, nnull, qraux, jpvt, q, ldq, nlaht,
&      r, ldr, nr, cr, ldcr, dr, lddr, wk, info)
      integer lds, nobs, nnull, jpvt(*), ldq, ldr, nr, ldcr, lddr, info
      double precision s(lds,*), qraux(*), q(ldq,*), nlaht, r(ldr,*), 
&     cr(ldcr,*), dr(lddr,*), wk(2,*)
      double precision dum, ddot
      integer i, j, n, n0
      info = 0
      if(.not.( nnull .lt. 0 .or. nnull .ge. nobs .or. nobs .gt. lds 
&     .or. nobs .gt. ldq .or. ldr .lt. nobs .or. nr .lt. 1 .or. ldcr 
&     .lt. nobs .or. lddr .lt. nnull ))goto 23000
      info = -1
      return
23000 continue
      n0 = nnull
      n = nobs - nnull
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, wk, 1)
      j=1
23002 if(.not.(j.le.nr))goto 23004
      if(.not.( n0 .gt. 0 ))goto 23005
      call dqrsl (s, lds, nobs, nnull, qraux, r(1,j), dum, r(1,j), dum, 
&     dum, dum, 01000, info)
23005 continue
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, wk, r(n0+2,j), dum, r(n0+
&     2,j), dum, dum, dum, 01000, info)
      j=j+1
      goto 23002
23004 continue
      call dset (n, 10.d0 ** nlaht, wk(2,1), 2)
      call daxpy (n, 1.d0, q(n0+1,n0+1), ldq+1, wk(2,1), 2)
      call dcopy (n-1, q(n0+1,n0+2), ldq+1, wk(1,2), 2)
      call dpbfa (wk, 2, n, 1, info)
      if(.not.( info .ne. 0 ))goto 23007
      info = -2
      return
23007 continue
      j=1
23009 if(.not.(j.le.nr))goto 23011
      call dpbsl (wk, 2, n, 1, r(n0+1,j))
      j=j+1
      goto 23009
23011 continue
      call dcopy (n-2, q(n0+2,n0+1), ldq+1, wk, 1)
      j=1
23012 if(.not.(j.le.nr))goto 23014
      call dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, wk, r(n0+2,j), r(n0+2,j),
&      dum, dum, dum, dum, 10000, info)
      j=j+1
      goto 23012
23014 continue
      j=1
23015 if(.not.(j.le.nr))goto 23017
      call dset (n0, 0.d0, cr(1,j), 1)
      call dcopy (n, r(n0+1,j), 1, cr(n0+1,j), 1)
      if(.not.( n0 .gt. 0 ))goto 23018
      call dqrsl (s, lds, nobs, nnull, qraux, cr(1,j), cr(1,j), dum, 
&     dum, dum, dum, 10000, info)
23018 continue
      j=j+1
      goto 23015
23017 continue
      if(.not.( n0 .gt. 0 ))goto 23020
      j=1
23022 if(.not.(j.le.nr))goto 23024
      i=1
23025 if(.not.(i.le.n0))goto 23027
      dr(i,j) = r(i,j) - ddot (n, r(n0+1,j), 1, q(n0+1,i), 1)
      i=i+1
      goto 23025
23027 continue
      call dtrsl (s, lds, n0, dr(1,j), 01, info)
      call dprmut (dr(1,j), n0, jpvt, 1)
      j=j+1
      goto 23022
23024 continue
23020 continue
      return
      end

