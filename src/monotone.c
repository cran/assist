#include <stdio.h>
#include <math.h>

/**   cubic on [0, T]  **/
static double
rc(double x, double y, long ord) {
  double tmp;
  tmp=(x+y-fabs(x-y))/2.0;
  if(ord==1)  tmp=(tmp)*(tmp)*(3.0*(x+y-tmp)-tmp)/6.0;
  return(tmp);
}


/** end of kernels  **/


void mono_rk(double *x, double *y, double *f,
        long *n1, long *n2, long *type, double *res)
{
  long i,j, t, s;
  double x1, y1, sum=0.0, sum_tmp;

  for(i=0;i< *n1; i++){
    x1=x[i+1]-x[i];
    sum=0.0;
    for(j=0; j< *n2; j++){
      y1=y[j+1]-y[j];
     sum_tmp= 0.2777778*0.2777778*(f[3*i]*f[3*j])*     rc(x[i]+x1*0.1127017, y[j]+y1*0.1127017, *type)+
              0.2777778*0.4444444*((f[3*i]*f[3*j+1])*  rc(x[i]+x1*0.1127017, y[j]+y1*0.5, *type      )+
 	           	           (f[3*i+1]*f[3*j])*  rc(x[i]+x1*0.5,       y[j]+y1*0.1127017, *type));
     sum_tmp+=0.4444444*0.4444444*((f[3*i+1]*f[3*j+1])*rc(x[i]+x1*0.5,       y[j]+y1*0.5, *type))+
              0.2777778*0.2777778*((f[3*i+2]*f[3*j+2])*rc(x[i]+x1*0.8872983, y[j]+ y1*0.8872983, *type));
     sum_tmp+=  0.2777778*0.2777778*((f[3*i]*f[3*j+2])*rc(x[i]+x1*0.1127017, y[j]+y1*0.8872983, *type)+
                                     (f[3*i+2]*f[3*j])*rc(x[i]+x1*0.8872983, y[j]+y1*0.1127017, *type))+
              0.4444444*0.2777778*((f[3*i+1]*f[3*j+2])*rc(x[i]+x1*0.5,       y[j]+y1*0.8872983, *type)+
                                   (f[3*i+2]*f[3*j+1])*rc(x[i]+x1*0.8872983, y[j]+y1*0.5, *type));  

      sum+=sum_tmp*x1*y1;  
      res[i*(*n2)+j]=sum;  
    }
  }
}

void mono_f(double *x, double *y, double *f,
        long *nx, long *ny, long *type, double *res)
{
  
  long i,j;
  double x1, y1, sum=0.0;

  for(i=0;i< *ny; i++){
    sum=0.0;    
    for(j=0; j< *nx; j++){
      x1=x[j+1]-x[j];
      sum+=x1*(0.2777778*(f[3*j]*  rc(x[j]+x1*0.1127017, y[i], *type)+
                      f[3*j+2]*rc(x[j]+x1*0.8872983, y[i], *type))+
           0.4444444* f[3*j+1]*rc(x[j]+x1*0.5,       y[i], *type));      
      res[i*(*nx)+j]=sum;  
    }
  }
}


void mono_s(double *f, double *x, long *n, double *res)
{
  long i;
  double sum=0.0;

  for(i=0;i< *n; i++){ 
    sum += (x[i+1]-x[i])*(0.2777778*(f[3*i]+f[3*i+2])+0.4444444*f[3*i+1]);      
    res[i]=sum;
  }
}
