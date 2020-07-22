#include "pluto.h"
#include "user_globals.h"
/* ***************************************************************** */
double LamCool (double T, double ZbyZsun)
/*!
 *   returns Lambda[T,Z] \equiv cooling rate per unit volume/nH^2 in cgs units.
 *   using Wiersma+ 2009 cooling table for CIE; not bothering about n_e/n_H 
 *   changing at lower temperatures for now as it is a small effect.
 *   imput parameters: T temperature in K; Z metallicity in solar units
 ******************************************************************* */
{
  int    klo, khi, kmid;
  static int ntab;
  double  Tmid, scrh, dT, Lambda_lo, Lambda_hi;
  static double *L_HHe_tab, *L_Z_tab, *T_tab;
  
  FILE *fcool;

/* -------------------------------------------
        Read tabulated cooling function
   ------------------------------------------- */

  if (T_tab == NULL){
    print (" > Reading table from disk...\n");
    fcool = fopen("CT_WSS09.dat","r");
    if (fcool == NULL){
      print ("! LamCool:  CT_WSS09.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_HHe_tab = ARRAY_1D(20000, double);
    L_Z_tab = ARRAY_1D(20000, double);
    T_tab = ARRAY_1D(20000, double);
    

    ntab = 0;
    while (fscanf(fcool, "%lf  %lf  %lf\n", T_tab + ntab, L_HHe_tab + ntab, L_Z_tab+ntab)!=EOF) {
      ntab++;
    }
  }

/* ---------------------------------------------
            Make sure that T is well-defined 
   --------------------------------------------- */

  if (T != T){
    printf (" ! Nan found in lam_cool \n");
    printf (" ! T = %12.6e\n", T);
    QUIT_PLUTO(1);
  }

  if (T < g_minCoolingTemp) { 
    return 0.0;
  }

/* ----------------------------------------------
        Table lookup by binary search  
   ---------------------------------------------- */

  klo = 0;
  khi = ntab - 1;

  if (T > T_tab[khi] || T < T_tab[klo]){
    //print (" ! T out of range   %12.6e\n",T);
    //QUIT_PLUTO(1);
    if (T>T_tab[khi]){
    	scrh = L_HHe_tab[khi]+ZbyZsun*L_Z_tab[khi];
	return scrh;
    }
    else{
	scrh = L_HHe_tab[klo]+ZbyZsun*L_Z_tab[klo];
	return scrh;	
    }
  }

  while (klo != (khi - 1)){
    kmid = (klo + khi)/2;
    Tmid = T_tab[kmid];
    if (T <= Tmid){
      khi = kmid;
    }else if (T > Tmid){
      klo = kmid;
    }
  }

/* -----------------------------------------------
    Compute and return Lambda 
   ----------------------------------------------- */

  dT       = T_tab[khi] - T_tab[klo];
  Lambda_lo   = L_HHe_tab[klo]+ZbyZsun*L_Z_tab[klo]; Lambda_hi = L_HHe_tab[khi]+ZbyZsun*L_Z_tab[khi];
  scrh     = Lambda_lo*(T_tab[khi] - T)/dT + Lambda_hi*(T - T_tab[klo])/dT;
  return scrh;
}

double MMWt_mu (double ZbyZsun,double XbyXsun)
{
    double Xsol = 0.7065, Ysol = 0.2806, Zsol = 1.-Xsol-Ysol, X, Y, Z;
    Z = Zsol*ZbyZsun;
    X = Xsol*XbyXsun;
    Y = 1.-X-Z;
    return 1./(2.*X+0.75*Y+9.*Z/16.);
}

double MMWt_muH (double ZbyZsun,double XbyXsun)
{
    double Xsol = 0.7065, Ysol = 0.2806, Zsol = 1.-Xsol-Ysol, X, Z;
    Z = Zsol*ZbyZsun;
    X = Xsol*XbyXsun; //(1.-Z)*Xsol/(1.-Zsol);
    return 1./X;
}

double MMWt_mui (double ZbyZsun,double XbyXsun)
{
    double Xsol = 0.7065, Ysol = 0.2806, Zsol = 1.-Xsol-Ysol, X, Y, Z;
    Z = Zsol*ZbyZsun;
    X = Xsol*XbyXsun;
    Y = 1.-X-Z;
    return 1./(X+Y/4.+Z/16.);
}

double MMWt_mue (double ZbyZsun,double XbyXsun)
{
    double Xsol = 0.7065, Ysol = 0.2806, Zsol = 1.-Xsol-Ysol, X, Y, Z;
    Z = Zsol*ZbyZsun;
    X = Xsol*XbyXsun;
    Y = 1.-X-Z;
    return 1./(X+0.5*Y+0.5*Z);
}
