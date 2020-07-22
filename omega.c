#include "pluto.h"
#include "lam_cool.h"
#include "background.h"
#include "user_globals.h"

double lam_T(double x)
/*!
 *   Calculate logarithmic derivative of cooling function wrt Temperature
 ******************************************************************* */
{
  
  
  double z_bg,x_bg,rho_bg,p_bg, T_bg, n_Hbg;
  double mue_bg=1.17, mu_bg=0.62;
  double gamma = 5./3.;
  double unit_q;
  
  double lam0,lam0p,dT, dlam0;
  
  rho_bg = rhobg(x);
  p_bg = pbg(x);
  
  unit_q = UNIT_DENSITY*pow(UNIT_VELOCITY,3.0); 
  unit_q /= UNIT_LENGTH;
  
  z_bg = zbg(x); 
  x_bg = xbg(x); 
  
  T_bg = p_bg/rho_bg*KELVIN*MMWt_mu(z_bg,x_bg);
  n_Hbg = rho_bg*UNIT_DENSITY/(MMWt_muH(z_bg,x_bg)*CONST_amu); 

  lam0 = LamCool(T_bg, z_bg);
  lam0 = lam0/unit_q;
  
  dT = 0.15*T_bg;
  
  lam0p = LamCool(T_bg+dT, z_bg);
  lam0p = lam0p/unit_q;
  
  dlam0 = lam0p-lam0;
  //printf("lam_T: %0.40lf\n",dlam0);

  return (T_bg/lam0)*(dlam0/dT);
  
  
}

double lam_z(double x)
/*!
 *   Calculate logarithmic derivative of cooling function wrt metallicity
 ******************************************************************* */
{
  
  
  double z_bg,x_bg,rho_bg,p_bg, T_bg, n_Hbg;
  double mue_bg=1.17, mu_bg=0.62;
  double gamma = 5./3.;
  double unit_q;
  
  double lam0,lam0p,dz, dlam0;
  
  rho_bg = rhobg(x);
  p_bg = pbg(x);
  
  unit_q = UNIT_DENSITY*pow(UNIT_VELOCITY,3.0); 
  unit_q /= UNIT_LENGTH;
  
  z_bg = zbg(x);
  x_bg = xbg(x); 
  
  T_bg = p_bg/rho_bg*KELVIN*MMWt_mu(z_bg,x_bg);
  n_Hbg = rho_bg*UNIT_DENSITY/(MMWt_muH(z_bg,x_bg)*CONST_amu); 

  lam0 = LamCool(T_bg, z_bg);
  lam0 = lam0/unit_q;
  
  dz = 0.001;
  
  lam0p = LamCool(T_bg, z_bg+dz);
  lam0p = lam0p/unit_q;
  
  dlam0 = lam0p-lam0;
  
  //printf("lam_z: %0.40lf\n",dlam0);
  
  return (z_bg/lam0)*(dlam0/dz);
  
  
}


/* ***************************************************************** */
double wi(double k, double x, int i)
/*!
 *   Calculate imaginary part of angular frequency
 *   i>0 : Isobaric
 *   i<=0 : Isochoric
 ******************************************************************* */
{
  
  double z_bg,x_bg,rho_bg,p_bg, T_bg, n_Hbg;
  double mue_bg=1.17, mu_bg=0.62;
  double gamma = 5./3.;
  double unit_q;
  
  double lam0, tc0, tti;
  
  rho_bg = rhobg(x);
  p_bg = pbg(x);
  
  unit_q = UNIT_DENSITY*pow(UNIT_VELOCITY,3.0); 
  unit_q /= UNIT_LENGTH;
  
  z_bg = zbg(x); 
  x_bg = xbg(x); 
  
  T_bg = p_bg/rho_bg*KELVIN*MMWt_mu(z_bg,x_bg);
  n_Hbg = rho_bg*UNIT_DENSITY/(MMWt_muH(z_bg,x_bg)*CONST_amu); 

  lam0 = LamCool(T_bg, z_bg);
  lam0 = n_Hbg*n_Hbg*lam0/unit_q;
  
  tc0 =  p_bg/(lam0*(gamma-1));
  
  if(i>0){
    tti = gamma*tc0/(2-lam_T(x));
    return 1/tti; 
  }
  else{
    return -1.0*lam_T(x)/tc0;
  }
}

/* ***************************************************************** */
double wr(double k, double x, int i)
/*!
 *   Calculate real part of angular frequency
 *   i>0 : Isobaric
 *   i<0 : Isochoric
 ******************************************************************* */
{
  if(i>0){
    double z_bg,x_bg,rho_bg,p_bg, T_bg,n_Hbg;
    double mu_ibg, mu_ebg, mu0_bg;
    double mue_bg=1.17, mu_bg=0.62;
    double gamma = 5./3.;
    double unit_q;
  
      double lam0, tc0;
      
      rho_bg = rhobg(x);
      p_bg = pbg(x);
      
      unit_q = UNIT_DENSITY*pow(UNIT_VELOCITY,3.0); 
      unit_q /= UNIT_LENGTH;
      
      z_bg = zbg(x);
      x_bg = xbg(x); 
      
      T_bg = p_bg/rho_bg*KELVIN*MMWt_mu(z_bg,x_bg);
      n_Hbg = rho_bg*UNIT_DENSITY/(MMWt_muH(z_bg,x_bg)*CONST_amu);
      
      mu_ibg = MMWt_mui(z_bg,x_bg);
      mu_ebg = MMWt_mue(z_bg,x_bg); 
      
      mu0_bg = 1/mu_ibg + 1/mu_ebg;
      mu0_bg = 1/mu0_bg;
      
      //printf("mu0: %lf|%lf\n--------------------------------------------------------         \n",mu0_bg,MMWt_muH(z_bg,x_bg));
    
      lam0 = LamCool(T_bg, z_bg);
      lam0 = n_Hbg*n_Hbg*lam0/unit_q;
      
      tc0 =  p_bg/(lam0*(gamma-1));
      
      //double m = m_z(x);
      
      double alpha,beta,epsilon,eta;
      
      alpha = (lam_z(x)/z_bg) - 3.*mu_ibg/16.;
      beta = 3.*mu_ibg/4. + mu_ebg/2.;
      
      epsilon = alpha + 3.*lam_T(x)*mu0_bg/16.;
      eta = beta + 5.*lam_T(x)*mu0_bg/4.;
      
      return (1./(gamma*k*tc0)) * (epsilon*m_z(x) + eta*m_x(x));
   }
   else{
      
      return 0.0;
   
   }
}

