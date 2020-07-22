#include "pluto.h"
#include "lam_cool.h"
#include "user_globals.h"

double rhobg(double x)
/*!
 *   Calculate logarithmic derivative of cooling function wrt Temperature
 ******************************************************************* */
{
  
  double rho_bg;
  double mue_bg=1.17, mu_bg=0.62;
  
  rho_bg = 0.1*(CONST_amu*mu_bg/UNIT_DENSITY);

  return rho_bg;
  
  
}

double pbg(double x)
/*!
 *   Calculate logarithmic derivative of cooling function wrt Temperature
 ******************************************************************* */
{
  double p_bg;
  
  
  p_bg = g_dp*CONST_kB*0.78*1.16e7/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);

  return p_bg;
  
  
}


double vbg(double x)
/*!
 *   Calculate logarithmic derivative of cooling function wrt Temperature
 ******************************************************************* */
{
  double v_bg;
  
  v_bg = 0.0;

  return v_bg;
  
  
}

double zbg(double x)
/*!
 *   Calculate logarithmic derivative of cooling function wrt Temperature
 ******************************************************************* */
{
  double m = g_m;   //defined in user_globals.h
  double c = g_cz;
  
  double z_bg;
  
  z_bg = m*x + c; 

  return z_bg;

}

double xbg(double x)
/*!
 *   Calculate logarithmic derivative of cooling function wrt Temperature
 ******************************************************************* */
{
  double m = g_mx;   //defined in user_globals.h
  double c = g_cx;
  
  double x_bg;
  
  x_bg = m*x + c; 

  return x_bg;

}

double m_z(double x)
/*!
 *   Calculate slope of metallicity distribution
 ******************************************************************* */
{
  double dx = 0.01;
  double dZ0 = zbg(x+dx) - zbg(x-dx);

  return 0.5*dZ0/dx;
}

double m_x(double x)
/*!
 *   Calculate slope of metallicity distribution
 ******************************************************************* */
{
  double dx = 0.01;
  double dX0 = xbg(x+dx) - xbg(x-dx);

  return 0.5*dX0/dx;
}

