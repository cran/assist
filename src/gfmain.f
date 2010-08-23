      subroutineqf(gf, alb,anc,n,irr,sigma,cc,lim1,acc,ith,trace,ifault)

      integerirr,lim1,ifault

      realsigma,cc,acc, gf

      realtrace(7),alb(irr),anc(irr)

      integern(irr),ith(irr)

      integerj,nj,nt,ntm

      realacc1,almx,xlim,xnt,xntm

      realutx,tausq,sd,aintv,aintv1,x,up,un,d1,d2,alj,ancj

      doubleprecisionaintl,ersm

      realpi,aln28,sigsq,almax,almin,amean,c

      integericount,ir,lim

      logicalndtsrt,fail

      common/qfcom/aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,icount,

     &ir,ndtsrt,fail,lim

      iround(x)=int(x+sign(0.5,x))

      pi=3.14159265358979

      aln28=.0866

      c=cc

      ir=irr

      lim=lim1

c      call intpr('a message', -1, ir, 1)

c      call dblepr or realpr 

      do 23000j=1,7

      trace(j)=0.0

23000 continue

23001 continue

      ifault=0

      icount=0

      aintl=0.0

      ersm=0.0

      gf=-1.0

      acc1=acc

      ndtsrt=.true.

      fail=.false.

      xlim=lim

      sigsq=sigma**2

      sd=sigsq

      almax=0.0

      almin=0.0

      amean=0.0

      continue

      j=1

23002 if(.not.(j.le.ir))goto 23004

      nj=n(j)

      alj=alb(j)

      ancj=anc(j)

      if(.not.(nj.lt.0.or.ancj.lt.0.0))goto 23005

      ifault=3

      goto100

23005 continue

      sd=sd+alj**2*(2*nj+4.0*ancj)

      amean=amean+alj*(nj+ancj)

      if(.not.(almax.lt.alj))goto 23007

      almax=alj

      goto 23008

23007 continue

      if(.not.(almin.gt.alj))goto 23009

      almin=alj

23009 continue

23008 continue

23003 j=j+1

      goto 23002

23004 continue

      if(.not.(sd.eq.0.0))goto 23011

      if(.not.(c.gt.0.0))goto 23013

      gf=1.0

      goto 23014

23013 continue

      gf=0.0

23014 continue

      goto100

23011 continue

      if(.not.(almin.eq.0.0.and.almax.eq.0.0.and.sigma.eq.0.0))goto 2301

     &5

      ifault=3

      goto100

23015 continue

      sd=sqrt(sd)

      if(.not.(almax.lt.-almin))goto 23017

      almx=-almin

      goto 23018

23017 continue

      almx=almax

23018 continue

      utx=16.0/sd

      up=4.5/sd

      un=-up

      callfindu(n,alb,anc,utx,.5*acc1)

      if(.not.(c.ne.0.0.and.almx.gt.0.07*sd))goto 23019

      tausq=.25*acc1/cfe(n,alb,anc,ith,c)

      if(.not.(fail))goto 23021

      fail=.false.

      goto 23022

23021 continue

      if(.not.(truncn(n,alb,anc,utx,tausq).lt..2*acc1))goto 23023

      sigsq=sigsq+tausq

      callfindu(n,alb,anc,utx,.25*acc1)

      trace(6)=sqrt(tausq)

23023 continue

23022 continue

23019 continue

      trace(5)=utx

      acc1=0.5*acc1

20    d1=ctff(n,alb,anc,acc1,up)-c

      if(.not.(d1.lt.0.0))goto 23025

      gf=1.0

      goto100

23025 continue

      d2=c-ctff(n,alb,anc,acc1,un)

      if(.not.(d2.lt.0.0))goto 23027

      gf=0.0

      goto100

23027 continue

      if(.not.(d1.gt.d2))goto 23029

      aintv=d1

      goto 23030

23029 continue

      aintv=d2

23030 continue

      aintv=2.0*pi/aintv

      xnt=utx/aintv

      xntm=3.0/sqrt(acc1)

      if(.not.(xnt.gt.xntm*1.5))goto 23031

      if(.not.(xntm.gt.xlim))goto 23033

      ifault=1

      goto100

23033 continue

      ntm=iround(xntm)

      aintv1=utx/xntm

      x=2.0*pi/aintv1

      if(.not.(x.le.abs(c)))goto 23035

      goto40

23035 continue

      tausq=cfe(n,alb,anc,ith,c-x)+cfe(n,alb,anc,ith,c+x)

      tausq=.33*acc1/(1.1*tausq)

      if(.not.(fail))goto 23037

      goto40

23037 continue

      acc1=.67*acc1

      callintegr(n,alb,anc,ntm,aintv1,tausq,.false.)

      xlim=xlim-xntm

      sigsq=sigsq+tausq

      trace(3)=trace(3)+1

      trace(2)=trace(2)+ntm+1

      callfindu(n,alb,anc,utx,.25*acc1)

      acc1=0.75*acc1

      goto20

23031 continue

40    trace(4)=aintv

      if(.not.(xnt.gt.xlim))goto 23039

      ifault=1

      goto100

23039 continue

      nt=iround(xnt)

      callintegr(n,alb,anc,nt,aintv,0.0,.true.)

      trace(3)=trace(3)+1

      trace(2)=trace(2)+nt+1

      gf=0.5-aintl

      trace(1)=ersm

      up=ersm

      x=up+acc/10.0

      continue

      j=1

23041 if(.not.(j.le.8))goto 23043

      if(.not.(j*x.eq.j*up))goto 23044

      ifault=2

23044 continue

23042 j=j*2

      goto 23041

23043 continue

100   trace(7)=icount

      return

      end

      subroutinecountr

      doubleprecisionaintl,ersm

      realpi,aln28,sigsq,almax,almin,amean,c

      integericount,ir,lim

      logicalndtsrt,fail

      common/qfcom/aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,icount,

     &ir,ndtsrt,fail,lim

      icount=icount+1

      if(.not.(icount.gt.lim))goto 23046

c      write(6,100)

c 100   format(' qf: cannot locate integration parameters'/)

      stop

23046 continue

      return

      end

      realfunctionalog1(x,first)

      realx

      logicalfirst

      reals,s1,term,y,ak

      f1(i)=s+term/ak

      if(.not.(abs(x).gt.0.1))goto 23048

      if(.not.(first))goto 23050

      alog1=alog(1.0+x)

      goto 23051

23050 continue

      alog1=alog(1.0+x)-x

23051 continue

      goto 23049

23048 continue

      y=x/(2.0+x)

      term=2.0*y**3

      ak=3.0

      if(.not.(first))goto 23052

      s=2.0

      goto 23053

23052 continue

      s=-x

23053 continue

      s=s*y

      y=y**2

      continue

      s1=f1(0)

23054 if(.not.(s1.ne.s))goto 23056

      ak=ak+2.0

      term=term*y

      s=s1

23055 s1=f1(0)

      goto 23054

23056 continue

      alog1=s

23049 continue

      return

      end

      realfunctionexp1(x)

      realx

      if(.not.(x.lt.-50.0))goto 23057

      exp1=0.0

      goto 23058

23057 continue

      exp1=exp(x)

23058 continue

      return

      end

      subroutineorder(alb,ith)

      realalb(ir)

      integerith(ir)

      integerj,k,k1,ithk

      realalj

      doubleprecisionaintl,ersm

      realpi,aln28,sigsq,almax,almin,amean,c

      integericount,ir,lim

      logicalndtsrt,fail

      common/qfcom/aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,icount,

     &ir,ndtsrt,fail,lim

      continue

      j=1

23059 if(.not.(j.le.ir))goto 23061

      alj=abs(alb(j))

      continue

      k=j-1

23062 if(.not.(k.gt.0))goto 23064

      ithk=ith(k)

      k1=k+1

      if(.not.(alj.gt.abs(alb(ithk))))goto 23065

      ith(k1)=ithk

      goto 23066

23065 continue

      goto20

23066 continue

23063 k=k-1

      goto 23062

23064 continue

      k=0

      k1=1

20    ith(k1)=j

23060 j=j+1

      goto 23059

23061 continue

      ndtsrt=.false.

      return

      end

      realfunctionerrbd(n,alb,anc,uu,cx)

      realu,uu,cx

      integern(ir)

      realalb(ir),anc(ir)

      realsum1,alj,ancj,x,y,const

      integerj,nj

      doubleprecisionaintl,ersm

      realpi,aln28,sigsq,almax,almin,amean,c

      integericount,ir,lim

      logicalndtsrt,fail

      common/qfcom/aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,icount,

     &ir,ndtsrt,fail,lim

      callcountr

      u=uu

      const=u*sigsq

      sum1=u*const

      u=2.0*u

      continue

      j=ir

23067 if(.not.(j.gt.0))goto 23069

      nj=n(j)

      alj=alb(j)

      ancj=anc(j)

      x=u*alj

      y=1.0-x

      const=const+alj*(ancj/y+nj)/y

      sum1=sum1+ancj*(x/y)**2

      sum1=sum1+nj*(x**2/y+alog1(-x,.false.))

23068 j=j-1

      goto 23067

23069 continue

      errbd=exp1(-0.5*sum1)

      cx=const

      return

      end

      realfunctionctff(n,alb,anc,accx,upn)

      realaccx,upn

      integern(ir)

      realalb(ir),anc(ir)

      realu1,u2,u,rb,const,c1,c2

      doubleprecisionaintl,ersm

      realpi,aln28,sigsq,almax,almin,amean,c

      integericount,ir,lim

      logicalndtsrt,fail

      common/qfcom/aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,icount,

     &ir,ndtsrt,fail,lim

      f1(i)=u2/(1.0+u2*rb)

      f2(i)=(c1-amean)/(c2-amean)

      u2=upn

      u1=0.0

      c1=amean

      if(.not.(u2.gt.0.0))goto 23070

      rb=almax

      goto 23071

23070 continue

      rb=almin

23071 continue

      rb=2.0*rb

      continue

      u=f1(0)

23072 if(.not.(errbd(n,alb,anc,u,c2).gt.accx))goto 23074

      u1=u2

      c1=c2

      u2=2.0*u2

23073 u=f1(0)

      goto 23072

23074 continue

      continue

      u=f2(0)

23075 if(.not.(u.lt.0.9))goto 23077

      u=(u1+u2)/2.0

      if(.not.(errbd(n,alb,anc,u/(1.0+u*rb),const).gt.accx))goto 23078

      u1=u

      c1=const

      goto 23079

23078 continue

      u2=u

      c2=const

23079 continue

23076 u=f2(0)

      goto 23075

23077 continue

      ctff=c2

      upn=u2

      return

      end

      realfunctiontruncn(n,alb,anc,uu,tausq)

      realu,uu,tausq

      integern(ir)

      realalb(ir),anc(ir)

      realsum1,sum2,prod1,prod2,prod3,alj,ancj,x,y,err1,err2

      integerj,nj,ns

      doubleprecisionaintl,ersm

      realpi,aln28,sigsq,almax,almin,amean,c

      integericount,ir,lim

      logicalndtsrt,fail

      common/qfcom/aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,icount,

     &ir,ndtsrt,fail,lim

      callcountr

      u=uu

      sum1=0.0

      prod2=0.0

      prod3=0.0

      ns=0

      sum2=(sigsq+tausq)*u**2

      prod1=2.0*sum2

      u=2.0*u

      continue

      j=1

23080 if(.not.(j.le.ir))goto 23082

      alj=alb(j)

      ancj=anc(j)

      nj=n(j)

      x=(u*alj)**2

      sum1=sum1+ancj*x/(1.0+x)

      if(.not.(x.gt.1.0))goto 23083

      prod2=prod2+nj*alog(x)

      prod3=prod3+nj*alog1(x,.true.)

      ns=ns+nj

      goto 23084

23083 continue

      prod1=prod1+nj*alog1(x,.true.)

23084 continue

23081 j=j+1

      goto 23080

23082 continue

      sum1=0.5*sum1

      prod2=prod1+prod2

      prod3=prod1+prod3

      x=exp1(-sum1-0.25*prod2)/pi

      y=exp1(-sum1-0.25*prod3)/pi

      if(.not.(ns.eq.0))goto 23085

      err1=1.0

      goto 23086

23085 continue

      err1=x*2.0/ns

23086 continue

      if(.not.(prod3.gt.1.0))goto 23087

      err2=2.5*y

      goto 23088

23087 continue

      err2=1.0

23088 continue

      if(.not.(err2.lt.err1))goto 23089

      err1=err2

23089 continue

      x=0.5*sum2

      if(.not.(x.le.y))goto 23091

      err2=1.0

      goto 23092

23091 continue

      err2=y/x

23092 continue

      if(.not.(err1.lt.err2))goto 23093

      truncn=err1

      goto 23094

23093 continue

      truncn=err2

23094 continue

      return

      end

      subroutinefindu(n,alb,anc,utx,accx)

      realutx,accx

      integern(ir)

      realalb(ir),anc(ir)

      realu,ut

      realdivis(4)

      integeri

      doubleprecisionaintl,ersm

      realpi,aln28,sigsq,almax,almin,amean,c

      integericount,ir,lim

      logicalndtsrt,fail

      common/qfcom/aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,icount,

     &ir,ndtsrt,fail,lim

      datadivis/2.0,1.4,1.2,1.1/

      ut=utx

      u=ut/4.0

      if(.not.(truncn(n,alb,anc,u,0.0).gt.accx))goto 23095

      continue

      u=ut

23097 if(.not.(truncn(n,alb,anc,u,0.0).gt.accx))goto 23099

      ut=ut*4.0

23098 u=ut

      goto 23097

23099 continue

      goto 23096

23095 continue

      ut=u

      continue

      u=u/4.0

23100 if(.not.(truncn(n,alb,anc,u,0.0).le.accx))goto 23102

      ut=u

23101 u=u/4.0

      goto 23100

23102 continue

23096 continue

      do 23103i=1,4

      u=ut/divis(i)

      if(.not.(truncn(n,alb,anc,u,0.0).le.accx))goto 23105

      ut=u

23105 continue

23103 continue

23104 continue

      utx=ut

      return

      end

      subroutineintegr(n,alb,anc,nterm,aintrv,tausq,main)

      integernterm

      realaintrv,tausq

      logicalmain

      integern(ir)

      realalb(ir),anc(ir)

      realainpi,u,sum1,sum2,sum3,x,y,z

      integerk,j,nj

      doubleprecisionaintl,ersm

      realpi,aln28,sigsq,almax,almin,amean,c

      integericount,ir,lim

      logicalndtsrt,fail

      common/qfcom/aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,icount,

     &ir,ndtsrt,fail,lim

      ainpi=aintrv/pi

      continue

      k=nterm

23107 if(.not.(k.ge.0))goto 23109

      u=(k+0.5)*aintrv

      sum1=-2.0*u*c

      sum2=abs(sum1)

      sum3=-0.5*sigsq*u**2

      continue

      j=ir

23110 if(.not.(j.gt.0))goto 23112

      nj=n(j)

      x=2.0*alb(j)*u

      y=x**2

      sum3=sum3-0.25*nj*alog1(y,.true.)

      y=anc(j)*x/(1.0+y)

      z=nj*atan(x)+y

      sum1=sum1+z

      sum2=sum2+abs(z)

      sum3=sum3-0.5*x*y

23111 j=j-1

      goto 23110

23112 continue

      x=ainpi*exp1(sum3)/u

      if(.not.(.not.main))goto 23113

      x=x*(1.0-exp1(-0.5*tausq*u**2))

23113 continue

      sum1=sin(0.5*sum1)*x

      sum2=0.5*sum2*x

      aintl=aintl+sum1

      ersm=ersm+sum2

23108 k=k-1

      goto 23107

23109 continue

      return

      end

      realfunctioncfe(n,alb,anc,ith,x)

      realx

      integern(ir),ith(ir)

      realalb(ir),anc(ir)

      realaxl,axl1,axl2,sxl,sum1,alj

      integerj,k,it,itk

      doubleprecisionaintl,ersm

      realpi,aln28,sigsq,almax,almin,amean,c

      integericount,ir,lim

      logicalndtsrt,fail

      common/qfcom/aintl,ersm,pi,aln28,sigsq,almax,almin,amean,c,icount,

     &ir,ndtsrt,fail,lim

      callcountr

      if(.not.(ndtsrt))goto 23115

      callorder(alb,ith)

23115 continue

      axl=abs(x)

      sxl=sign(1.0,x)

      sum1=0.0

      continue

      j=ir

23117 if(.not.(j.gt.0))goto 23119

      it=ith(j)

      if(.not.(alb(it)*sxl.gt.0.0))goto 23120

      alj=abs(alb(it))

      axl1=axl-alj*(n(it)+anc(it))

      axl2=alj/aln28

      if(.not.(axl1.gt.axl2))goto 23122

      axl=axl1

      goto 23123

23122 continue

      if(.not.(axl.gt.axl2))goto 23124

      axl=axl2

23124 continue

      sum1=(axl-axl1)/alj

      continue

      k=j-1

23126 if(.not.(k.gt.0))goto 23128

      itk=ith(k)

      sum1=sum1+(n(itk)+anc(itk))

23127 k=k-1

      goto 23126

23128 continue

      goto10

23123 continue

23120 continue

23118 j=j-1

      goto 23117

23119 continue

10    if(.not.(sum1.gt.100.0))goto 23129

      cfe=1.0

      fail=.true.

      goto 23130

23129 continue

      cfe=exp(log(2.0)*(sum1/4.0))/(pi*axl**2)

23130 continue

      return

      end





















