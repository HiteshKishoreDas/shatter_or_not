#include "pluto.h"
#include "lam_cool.h"
#include "background.h"
/* ***************************************************************** */
void LamHeatxt (Grid *GXYZ)
/*!
 *   creates txt file for heating rate per unit volume/nH^2 in cgs units.
 *   Calculated using Wiersma+ 2009 cooling table for CIE from background 
 *   values; not bothering about n_e/n_H changing at lower temperatures 
 *   for now as it is a small effect.
 *   imput parameters: Grid Structure, index of position
 ******************************************************************* */
{
  
  static int ntab;
  
  static double *nH_tab, *heat_tab;
  static int *i_tab;
  double unit_q;
  
  double ini_heat;                // heating source
  double T_bg;                    // background temp
  double *x1;                     // x coordinate
  double mu_bg =0.62,n_Hbg;
  
  FILE *fheatr,*fheatw;
  
  x1 = GXYZ->x[IDIR];
  unit_q = UNIT_DENSITY*pow(UNIT_VELOCITY,3.0); 
  
  unit_q /= UNIT_LENGTH;
 
  int filexist = 0;
  
  if((fheatr = fopen("heatfunc.txt","r")))
  {
    fclose(fheatr);
    filexist = 1;
  }
  
  if (!filexist)
  {
    print ("! LamHeat:  Creating heatfunc.txt ....\n");
    fheatw = fopen("heatfunc.txt","w+");
    int i;

    for (i = 0;i<NX1;i++)
    {    
        T_bg = pbg(x1[i+2])/rhobg(x1[i+2])*KELVIN*MMWt_mu(zbg(x1[i+2]),xbg(x1[i+2]));
        n_Hbg = rhobg(x1[i+2])*UNIT_DENSITY/(MMWt_muH(zbg(x1[i+2]),xbg(x1[i+2]))*CONST_amu); 
    
        ini_heat= LamCool(T_bg, zbg(x1[i+2]));
        ini_heat = n_Hbg*n_Hbg*ini_heat/unit_q;
        printf("%lf %d|\n",ini_heat,i); 

        fprintf(fheatw, "%0.50lf  %d  %0.50lf\n",n_Hbg,i, ini_heat);
        
    }
    fclose(fheatw);

    printf("LamHeatxt(): heatfunc.txt created......\n");
  }
   
}
