#include "pluto.h"
#include "lam_cool.h"
#include "omega.h"
#include "background.h"
#include "user_globals.h"


//#####################################################################
// Velocity perturbation amplitude
//#####################################################################


double delvr(double k,double x, double del_rho)
/*!
 *   Calculate real part of velocity perturbation amplitude
 ******************************************************************* */
{
  double del_vr;
  int iso = g_iso;
  
  del_vr = (wr(k,x,iso)/(k*rhobg(x))) * del_rho; 

  return del_vr; 
}

double delvi(double k,double x, double del_rho)
/*!
 *   Calculate imaginary part of velocity perturbation amplitude
 ******************************************************************* */
{
  double del_vi;
  int iso = g_iso;
  
  del_vi = (wi(k,x,iso)/(k*rhobg(x))) * del_rho; 

  return del_vi;  
}


//#####################################################################
// Pressure perturbation amplitude
//#####################################################################


double delpr(double k,double x, double del_rho)
/*!
 *   Calculate real part of pressure perturbation amplitude
 ******************************************************************* */
{
  double del_pr;
  int iso = g_iso;
  
  del_pr = ( ( pow(wr(k,x,iso),2.0) - pow(wi(k,x,iso),2.0) )/ (k*k) )* del_rho; 

  return del_pr; 
}

double delpi(double k,double x, double del_rho)
/*!
 *   Calculate imaginary part of pressure perturbation amplitude
 ******************************************************************* */
{
  double del_pi;
  int iso = g_iso;
  
  del_pi = ( 2.*wr(k,x,iso)*wi(k,x,iso) / (k*k) )* del_rho; 

  return del_pi; 
}


//#####################################################################
// Metallicity perturbation amplitude
//#####################################################################


double delzi(double k,double x, double del_rho)
/*!
 *   Calculate imaginary part of metallicity perturbation amplitude
 ******************************************************************* */
{
  double del_zi;  
  double mz = m_z(x);
  
  del_zi = ( -1.*mz/ ( k * rhobg(x) ) ) * del_rho; 

  return del_zi; 
}

//#####################################################################
// H mass fraction perturbation amplitude
//#####################################################################


double delxi(double k,double x, double del_rho)
/*!
 *   Calculate imaginary part of X perturbation amplitude
 ******************************************************************* */
{
  double del_xi;  
  double mx = m_x(x);
  
  del_xi = ( -1.*mx/ ( k * rhobg(x) ) ) * del_rho; 

  return del_xi; 
}


//#####################################################################
// Magnitude and argument of a complex number
//#####################################################################


double mag(double a, double b)
/*!
 *   Calculate magnitude of complex amplitude
 ******************************************************************* */
{
    double mag_c;
    
    mag_c = sqrt( a*a + b*b );
    
    return mag_c;
}

double arg(double a, double b)
{
    double arg_c;
    
    if (a == 0.0)
    {
        if (b>0.0) {return CONST_PI/2.;}
        else if (b<0.0) {return -1.*CONST_PI/2.;}
        else 
        {
            return 0.0;
            //printf("Im() = %lf\n",b);
            //printf("! eigenmode.c:arg()-Invalid input...\n");
            //QUIT_PLUTO(1);
        }
    }
    
    arg_c = atan(b/a);
    
    return arg_c;
} 



