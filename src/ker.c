#include <math.h>
#include <stdio.h>

/***~~~  cubic on [0, 1] ~~~***/
static double
dk4(double x) {
  double val;
  x=fabs(x); 
  val = x - 0.5;
  val *= val;
  val = (val * val - val/2. + 7./240.)/24.;
  return(val);
}

static double
dk2(double x) {
   double value;
   x= fabs(x);
   value= x - 0.5 ;
   value *= value;
   value= (value-1./12.)/2.;
   return(value);
}

static double                              
rc(double x, double y) {                   
   double value;                              
   value = dk2(x) * dk2(y)- dk4 (x - y);     
   return(value);                               
   }      

void
cubic_ker1(double *x, double *y, long *N, long *M, double *val){            
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = rc(x[i], y[j]);
    }
  }
}

/*~~~~   cubic on [0, T] ~~~~*/
static double
rkc(double x, double y) {
  double val, tmp;
  tmp=(x+y-fabs(x-y))/2.0;
  val=(tmp)*(tmp)*(3.0*(x+y-tmp)-tmp)/6.0;
  return(val);
}

void
cubic_ker2(double *x, double *y, long *N, long *M, double *val){
  /* this function is to calculate the rk for cubic spline  */
  /* with domain [0, T].                                    */
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = rkc(x[i], y[j]);
    }
  }
}


/***~~~~  periodic (m=1) ~~~~***/
void
period_ker(double *x, double *y, long *N, long *M, double *val){

  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = -dk4(x[i]-y[j]);
    }
  }
}

/***~~~~ period (m-order) ***/
static double
per_rkm(double x, long m) {
  double val=0;
  long i;

  for(i=2; i<=50; i++)
     val += (2.0/(pow(cos(6.283 * i), 2*m)))*cos(6.282*x);

  return(val);
}

void
mperiod_ker(double *x, double *y, long *N, long *M, long *ord, double *val){
  /* this function is to calculate the rk for cubic spline  */
  /* with domain [0, T].                                    */
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = per_rkm(x[i]-y[j], *ord);
    }
  }
}


/***~~~   expLspline  ~~~***/
static double
rk_exp(double x, double y, double r) {
  double tmp, val;
  
  tmp=(x+y- fabs(x-y))/2.0;
  val = r*tmp+exp(-r*x)+exp(-r*y)-0.5*exp(-r*x-r*y);
  val +=0.5*exp(r*(2.0*tmp-x-y))-exp(r*(tmp-y))-exp(r*(tmp-x));
  val= val/r/r/r;
  return(val);
 }

void
expLspline_ker(double *x, double *y, double *r, long *N, long *M, double *val) {
  long  i, j;

  for(i=0; i< *N; i++) 
      for(j=0; j< *M; j++) val[i*(*M)+j]=rk_exp(x[i], y[j], *r);
}

/***~~~~    sinLspline ~~~~***/
static double
rk_sinL0(double x) {
  double val=0;
  long i;

  for(i=2; i<=50; i++)
     val += 2 * cos(6.283 * i * x)/1558.545/(i*i*2-1)/(i*i*2-1);

  return(val);
}

void
sinLspline_ker0(double *x, double *y, long *N, long *M, double *val){
  /* this function is to calculate the rk for cubic spline  */
  /* with domain [0, T].                                    */
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = rk_sinL0(x[i]-y[j]);
    }
  }
}
static double
rk_sinL1(double x) {
  double val=0;
  long i;

  for(i=2; i<=50; i++)
     val += 2 * cos(6.283 * i * x)/61528.9/(i*i*2-1)/(i*i*2-1)/i/i;

  return(val);
}

void
sinLspline_ker1(double *x, double *y, long *N, long *M, double *val){
  /* this function is to calculate the rk for cubic spline  */
  /* with domain [0, T].                                    */
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = rk_sinL1(x[i]-y[j]);
    }
  }
}


static double
rk_sinL4p(double x) {
  double val=0;
  long i;

  for(i=3; i<=50; i++)
      /* val += 2.0 * cos(6.283 * i * x)/(i*i*2-1)/(i*i*2-1)/(i*i-4.0)/(i*i-4.0);*/
     val += 2.0 * cos(6.283 * i * x)/256.0/9488.531/(i*i*2-1)/(i*i*2-1)/(i*i-4.0)/(i*i-4.0);

  return(val);
}

void
sinLspline_ker4p(double *x, double *y, long *N, long *M, double *val){
  /* this function is to calculate the rk for cubic spline  */
  /* with domain [0, T].                                    */
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = rk_sinL4p(x[i]-y[j]);
    }
  }
}


/***~~~~  factor ~~~~***/
void
factor_ker(long *x, long *y, long *N, long *M, long *val){

  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = (x[i]==y[j]);
    }
  }
}

/***~~~~ linear-periodic ~~~~***/
static double
rk_linP(double x, double y){
   double val, tmp;
   tmp=(x+y-fabs(x-y))/2.0;
   val= tmp*tmp*tmp/3.0-(x+y)*tmp*tmp/2.0 + x*y*tmp;
   val += x*cos(y)-sin(y)-cos(tmp-y)*(x-tmp) -sin(tmp-y);
   val += y*cos(x)-sin(x)-cos(tmp-x)*(y-tmp) -sin(tmp-x);
   val += tmp/2.0 * cos(y-x) -0.25*sin(2*tmp-x-y) - 0.5*sin(x+y);
   return(val);
 }

void
  linPeriod_ker(double *x, double *y, long *N, long *M, double *val){
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = rk_linP(x[i], y[j]);
    }
  }
}

/***~~~   logistic   ~~~~***/
static double
rk_logit(double x, double y) {
  double val, tmp;
  tmp=(x+y-fabs(x-y))/2.0;
  val= tmp-2.0*exp(-tmp)-0.5*exp(-2.0*tmp)+2.5;
  val = val *exp(x+y)/(1+exp(x))/(1+exp(y));
  return(val);
}
void
logit_ker(double *x, double *y, long *N, long *M, double *val){
  /* this function is to calculate the rk for cubic spline  */
  /* with domain [0, T].                                    */
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = rk_logit(x[i], y[j]);
    }
  }
}

/***~~~   quintic [0, T]  ~~~~***/
static double
rk_quint(double x, double y) {
  double val, tmp;
  tmp=(x+y-fabs(x-y))/2.0;
  val=(pow(x, 5)-pow(x-tmp, 5))/20.0+(y-x)/8.0 *(pow(x, 4)-pow(x-tmp, 4));
  val += (y-x)*(y-x)*(pow(x, 3)-pow(x-tmp, 3))/12.0;
  return(val);
}

void
quintic_ker2(double *x, double *y, long *N, long *M, double *val){
  /* this function is to calculate the rk for quint spline  */
  /* with domain [0, T].                                    */
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = rk_quint(x[i], y[j]);
    }
  }
}

/***=== quintic spline [0, 1] ===***/
static double
quint_k3(double x){
  double val;
  x=fabs(x); 
  val = x - 0.5;
  val= val*val*val-0.25*val;
  return(val/6.0);
}
static double
quint_k6(double x){
  double val;
  x=fabs(x); 
  val = x - 0.5;
  val *=val;
  val= val*val*val-1.25*val*val+7.0/16.0*val-31.0/1344.0;
  return(val/720.0);
}    

static double
quint_rk(double x, double y){
  double val;
  val=quint_k3(x)*quint_k3(y)+quint_k6(x-y);
  return(val);
}

void
quintic_ker1(double *x, double *y, long *N, long *M, double *val){
  /* this function is to calculate the rk for quint spline  */
  /* with domain [0, 1].                                    */
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = quint_rk(x[i], y[j]);
    }
  }
}


/***~~~   septic [0 T]  ~~~~***/
static double
rk_sept(double x, double y) {
  double val, tmp;
  tmp=(x+y-fabs(x-y))/2.0;
  val=(pow(x, 7)-pow(x-tmp, 7))/7.0+(y-x)/2.0 *(pow(x, 6)-pow(x-tmp, 6));
  val +=0.6*(y-x)*(y-x)*(pow(x, 5)-pow(x-tmp, 5))+0.25*pow(y-x, 3)*(pow(x, 4)-pow(x-tmp, 4));
  return(val);
}

void
septic_ker2(double *x, double *y, long *N, long *M, double *val){                          
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = rk_sept(x[i], y[j]);
    }
  }
}

/***=== septic spline [0, 1] ===***/
static double
sep_k4(double x){
  double val;
  x=fabs(x); 
  val = x - 0.5;
  val *=val;
  val= val*val-0.5*val+7.0/240.0;
  return(val/24.0);
}
static double
sep_k8(double x){
  double val;
  x=fabs(x); 
  val = x - 0.5;
  val *=val;
  val= val*val*val*val-7.0/3.0*val*val*val+49.0/24.0*val*val-31.0/48.0*val+127.0/3840.0;
  return(val/720.0/56.0);
}    

static double
sep_rk(double x, double y){
  double val;
  val=sep_k4(x)*sep_k4(y)-sep_k8(x-y);
  return(val);
}

void
septic_ker1(double *x, double *y, long *N, long *M, double *val){
  /* this function is to calculate the rk for septic spline  */
  /* with domain [0, 1].                                    */
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = sep_rk(x[i], y[j]);
    }
  }
}

/***~~~~   thin-plate ~~~***/
void
tp_ker(double *x, double *y, long *d, long *m, long *N, long *M , double *val){

  long i, j, k;
  double tmp;

 for(i=0; i< *N; i++) 
      for(j=0; j< *M; j++){
         tmp=0;
         for(k=0; k< *d; k++) 
            tmp+= (x[i+k*(*N)]-y[j+k*(*M)])*(x[i+k*(*N)]-y[j+k*(*M)]);         
         val[i*(*M)+j]= pow(tmp, (*m) - (*d)/2.0);
         if(((*d)%2)==0) {             
               if(tmp>0) val[i*(*M)+j]=  val[i*(*M)+j] * log(tmp);
         }
    }
}


void
tp_term(long *dm, long *order, long *p){

long i, j, k, sum;
long pos=0;

 for(i=0; i< pow( (*order), *dm); i++){
   k=i;
   sum=0;
   for(j=0; j< *dm; j++){
       p[pos+j]= k % (*order);
       sum += p[pos+j];
       k= k / (*order);
   }
   if(sum< (*order)) pos += (*dm);
 }
}


/***~~~~   linear    ~~~~*/
static double
dk2_lin(double x) {
   double value;
   x= fabs(x);
   value= x - 0.5 ;
   value *= value;
   value= (value-1./12.)/2.;
   return(value);
}

static double                              
rc_lin(double x, double y) {                   
   double value;                              
   value = (x-0.5) * (y-0.5) +dk2_lin(x - y);     
   return(value);                               
   }

void
linear_ker1(double *x, double *y, long *N, long *M, double *val){
  /* this function is to calculate the rk for linear spline  */
  /* with domain [0, 1].                                    */
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = rc_lin(x[i], y[j]);
    }
  }
}

/***  linear [0, T]  ***/
static double
rk_lin(double x, double y){
  double val;
  val=x;
  if(x>y) val=y;
  return(val);
}

void
linear_ker2(double *x, double *y, long *N, long *M, double *val){
  /* this function is to calculate the rk for linear spline  */
  /* with domain [0, 1].                                    */
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = rk_lin(x[i], y[j]);
    }
  }
}
/*** ==== following are derivatives of rks ====***/
/***===  periodic  ===***/
static double
dkp(double x) {
  long i;
  double val=0;

  for(i=1; i<=100; i++)
      val += sin(6.283 * i * x)/124.0251/i/i/i;
  return(val);
}

void
dperiod_ker(double *x, double *y, long *N, long *M, double *val){
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = dkp(x[i]-y[j]);
    }
  }
}

/***~~~~    cubic on [0,1] ~~~~***/
static double
dkc(double x, double y) {
  double val, dd;
 
  dd=fabs(x-y);  
  val= ((dd-0.5)*(dd-0.5)*(dd-0.5)-(dd-0.5)/4.0)/6.0;
  if(x<y) return(val);
  else return(-1.0*val);
}

void
dcubic_ker1(double *x, double *y, long *N, long *M, double *val){
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = dkc(x[i], y[j]);
    }
  }
}

/***===   cubic on [0, T]  ===***/
static double
dkc2(double x, double y) {
  double val;
  
  if(x<y) val= x*x*0.5;
  else val= ( 2.0*x*y-y*y )/2.0;
  return(val);
}

void
dcubic_ker2(double *x, double *y, long *N, long *M, double *val){
  long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = dkc2(x[i], y[j]);
    }
  }
}

/***===  expLspline  ===***/
static double
dkexp(double x, double y, double r) {
  double tmp, val, tmp2=0.0;

  if( y< x ) tmp2=1.0;
  tmp=(x+y- fabs(x-y))/2.0;
  val = tmp2-exp(-r*y)+0.5*exp(-r*x-r*y);
  val +=(0.5-tmp2)*exp(r*(2.0*tmp-x-y));
  val =val/r/r;
  return(val);
 }

void
dexpLspline_ker(double *x, double *y, double *r, long *N, long *M, double *val)
{
  long  i, j;
  for(i=0; i< *N; i++)
      for(j=0; j< *M; j++) val[i*(*M)+j]=dkexp(x[i], y[j], *r);
}   

/***===  sinLspline  ===***/
static double
dksin(double x) {
  double val=0;
  long i;

  for(i=2; i<=50; i++)
     val +=  0.008062884*i * sin(6.283 * i * x)/(i*i*2-1)/(i*i*2-1);

  return(val);
}

void
dsinLspline_ker(double *x, double *y, long *N, long *M, double *val){
   long i,j;
   for(i = 0; i < *N; i++) {
    for(j = 0; j < *M; j++) {
      val[(i*(*M))+j] = dksin(x[i]-y[j]);
    }
  }
}    

/*  Sphere Kernel */
double cos_angle(double a1, double b1, double a2, double b2)
{
  double val;
  /*  val=sin(a1)*cos(b1)*sin(a2)*cos(b2)+cos(a1)*cos(b1)*cos(a2)*cos(b2)+
      sin(b1)*sin(b2); */
  val=cos(b1)*cos(b2)*cos(a1-a2)+sin(b1)*sin(b2);
  return val;
}

double rk_q(double z, long ord){
  double w, a, c, val;
  if(z>0.9999999){
    switch(ord){
     case 2:
       val=0.5;
       break;
    case 3:
      val=0.25;
      break;
    case 4:
      val=1.0/6.0;
      break;
    case 5:
      val=105.0/840.0;
      break;
    case 6:
      val=0.1;
      break;
    default:
      val=-1;
    }
  }
  else{
     w=(1-z)/2.0;
     c=2*sqrt(w);
     a=log(1.0+1.0/sqrt(w));
     switch(ord){
      case 2:
          val=(a*(12.0*w*w-4.0*w)-6.0*c*w+6.0*w+1)/2.0;
          break;
      case 3:
         val=(a*(840*w*w*w*w-720*w*w*w+72*w*w)+420*w*w*w+
                  c*(220*w*w-420*w*w*w)-150*w*w-4*w+3)/12;
         break;
     case 4:
       val=(a*(27720.0*pow(w, 6.0)-37800.0*pow(w, 5.0)+12600*pow(w, 4.0)-600*w*w*w)+
	    13860*pow(w, 5.0)+c*(-13860*pow(w, 5.0)+14280.0*pow(w, 4.0)-2772*w*w*w)-
            11970*pow(w, 4.0)+1470*w*w*w+15*w*w-3*w+5.0)/30.0;
       break;
     case 5:
       val=(a*(10810800.0*pow(w, 8.0)-20180160.0*pow(w, 7.0)+11642400.0*pow(w, 6)-
	       2116800.0*pow(w, 5.0)+58800.0*pow(w, 4.0))+5405400.0*pow(w, 7.0)+
            c*(-5405400.0*pow(w,7.0)+8288280.0*pow(w, 6.0)-3538920.0*pow(w, 5.0)+
	       363816.0*pow(w, 4.0))-7387380.0*pow(w, 6.0)+2577960.0*pow(w, 5.0)-
	    159810*pow(w, 4.0)-840*w*w*w+84*w*w-40*w+105)/840.0;
       break;
     case 6:
       val=(a*(232792560.0*pow(w, 10.0)-551350800.0*pow(w, 9.0)+454053600.0*pow(w, 8.0)-
	       151351200.0*pow(w, 7.0)+17463600.0*pow(w, 6.0)-317520.0*pow(w, 5.0))+
	    116396280.0*pow(w, 9.0)+
           c*(-116396280*pow(w, 9.0)+236876640.0*pow(w, 8.0)-
	    158414256.0*pow(w, 7.0)+38507040.0*pow(w, 6.0)-2462680.0*pow(w, 5.0))-
	    217477260.0*pow(w, 8.0)+127987860.0*pow(w, 7.0)-24954930.0*pow(w, 6.0)+
            930006.0*pow(w, 5.0)+2940.0*pow(w, 4.0)-180.0*w*w*w+45*w*w-35*w+126.0)/1260.0;
       break;
     default:
       val=-1;
     }    
  }
  return val;
 }
  
void sphere_ker(double *x1, double *y1, double *x2, double *y2, long *len1, 
               long *len2, long *ord, double *rk)
{
 long i, j;
 
 for(i=0; i< *len1; i++)
    for(j=0; j< *len2; j++){
      rk[j+i*(*len2)]=rk_q(cos_angle(x1[i], y1[i], x2[j], y2[j]), *ord)-1.0/(2*(*ord)-1);
     }
}

