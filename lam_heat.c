#include "pluto.h"

/* ***************************************************************** */
double LamHeat (int xi)
/*!
 *   returns equiv heating rate per unit volume/nH^2 in cgs units.
 *   Calculated using Wiersma+ 2009 cooling table for CIE from background 
 *   values; not bothering about n_e/n_H changing at lower temperatures 
 *   for now as it is a small effect.
 *   imput parameters: Grid Structure, index of position
 ******************************************************************* */
{
  static int ntab;
  static double *nH_tab, *heat_tab;
  static int *i_tab;
  
  FILE *fheatr;
  
  

  if (heat_tab == NULL)
  {
    print (" > Reading table from disk...\n");
    fheatr = fopen("heatfunc.txt","r");
    nH_tab = ARRAY_1D(20000, double);
    heat_tab = ARRAY_1D(20000, double);
    i_tab = ARRAY_1D(20000, int);

    
    ntab = 0;
    while (fscanf(fheatr, "%lf  %d  %lf\n", nH_tab + ntab, i_tab + ntab, heat_tab + ntab)!=EOF)
    {
      ntab++;
    }
    fclose(fheatr);
  }
  
  if (xi>=NX1 || xi < 0) {printf("! LamHeat:  xi out of range");}
  
  return heat_tab[xi]; 

}

