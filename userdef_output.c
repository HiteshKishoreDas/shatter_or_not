#include "pluto.h"
#include "omega.h"
#include "user_globals.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;  

  double ***shock, ***wreal, ***wimg, ***lmT,***lmZ;
  double kw = g_kw; //Defined in user_globals.h
  int iso = g_iso;
  double *x1;
  
  x1 = grid->x[IDIR];
  
  shock = GetUserVar("sk");
  wreal = GetUserVar("wr");
  wimg = GetUserVar("wi");
  lmT = GetUserVar("lamT");
  lmZ = GetUserVar("lamZ");

  DOM_LOOP(k,j,i){

	  shock[k][j][i] = d->flag[k][j][i];
	  wreal[k][j][i] = wr(kw,x1[i],iso);
	  wimg[k][j][i] = wi(kw,x1[i],iso);
	  lmT[k][j][i] = lam_T(x1[i]);
	  lmZ[k][j][i] = lam_z(x1[i]);
  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

#ifdef PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}





