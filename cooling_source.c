/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Simplify original PLUTO cooling routine to be consistent with the 
  semi-implicit treatment of Sharma et al. 2010.

  On output, a time-step estimate for the next time level is computed using
  the relative or absolute variation obtained during the integration of the ODE 
  (from t(n) --> t(n+1))  
  \f[
     \Delta t_c = \min_{ijk}\left[\frac{\Delta t^n M_R}{\epsilon}\right]
     \,\quad\rm{where}\quad
     \epsilon = \max\left(\left|\frac{p^{n+1}}{p^n} - 1\right|,\,
                |X^{n+1}-X^n|\right)
  \f]
  where \f$ M_R \f$ is the maximum cooling rate (defined by the global variable  
  ::g_maxCoolingRate) and X are the chemical species.
  
  \b References
     - "THERMAL INSTABILITY WITH ANISOTROPIC THERMAL CONDUCTION AND ADIABATIC 
     COSMIC RAYS: IMPLICATIONS FOR COLD FILAMENTS IN GALAXY CLUSTERS" \n 
     Sharma, Parrish and Quataert, ApJ (2010) 720, 652
     
  \authors A. Mignone (mignone@ph.unito.it)\n
           P. Sharma 
  \date    November 12, 2018
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "lam_heat.h"
#include "lam_heatxt.h"
#include "lam_cool.h"
#include "user_globals.h"

/* ********************************************************************* */
void CoolingSource (const Data *d, double dt, timeStep *Dts, Grid *GXYZ)
/*!
 * Integrate cooling and reaction source terms.
 *
 * \param [in,out]  d   pointer to Data structure
 * \param [in]     dt   the time step to be taken
 * \param [out]    Dts  pointer to the timeStep structure
 * \param [in]    GXYZ  pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int k, j, i, ncool, cnt;
  double scrh, dtsub;
  double T0, n_H;
  double q_m, q_p, vol, dq, unit_q;
  double mu, mue, mui;
  double sendarray[2], recvarray[2];
  double LamCool(double,double), MMWt_mu(double,double), MMWt_muH(double,double);
  double xloc, xglob;
  double *x1;
  double noch = g_noch;
  double xend,xbeg;
  double ln;
/*  ----------------------------------------------------------- 
                   first calculate the volume avg n_H^2 Lambda 
    -----------------------------------------------------------  */
  
  q_m = 0.0; vol = 0.0; unit_q = UNIT_DENSITY*pow(UNIT_VELOCITY,3.0); unit_q /= UNIT_LENGTH;

  x1= GXYZ->x[IDIR];
    
  xbeg = x1[IBEG];
  //ln = x1[1];
  xend = x1[IEND];
  
  //printf("%f-------------\n",xbeg);
  //printf("%f@@@@@@@@@@@@\n",ln);
  //printf("%f ##############\n",xend);
    
  LamHeatxt(GXYZ);

  DOM_LOOP(k,j,i){  /* -- span the computational domain -- */

    T0  = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*KELVIN*MMWt_mu(d->Vc[TRC][k][j][i],d->Vc[TRC+1][k][j][i]); //in cgs units
    scrh = LamCool(T0,d->Vc[TRC][k][j][i]);
    n_H = d->Vc[RHO][k][j][i]*UNIT_DENSITY/(MMWt_muH(d->Vc[TRC][k][j][i],d->Vc[TRC+1][k][j][i])*CONST_amu); 
    scrh = n_H*n_H*scrh/unit_q; //code units
  /* ------------------------------------------
      Suggest next time step
     ------------------------------------------ */
    scrh = d->Vc[PRS][k][j][i]/(scrh*(g_gamma-1.0)); //code units
    Dts->dt_cool = MIN(Dts->dt_cool, scrh);

  } /* -- end loop on points -- */

  #ifdef PARALLEL
  xloc = Dts->dt_cool;
  MPI_Allreduce (&xloc, &xglob, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  Dts->dt_cool = xglob;
  #endif

  ncool = ceil(dt/Dts->dt_cool); ncool = MAX(ncool, 50); // imposed maximum subcycling of 50
  dtsub = dt/ncool;
  
  
  for (cnt=0; cnt<ncool; ++cnt) {
    q_p = 0.0; vol = 0.0;
    DOM_LOOP(k,j,i){  /* -- span the computational domain -- */ //first apply cooling
      T0  = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*KELVIN*MMWt_mu(d->Vc[TRC][k][j][i],d->Vc[TRC+1][k][j][i]); //in cgs units
      //printf("%f \t %0.50f\n", T0,d->Vc[PRS][k][j][i]);
      scrh = LamCool(T0,d->Vc[TRC][k][j][i]);
      n_H = d->Vc[RHO][k][j][i]*UNIT_DENSITY/(MMWt_muH(d->Vc[TRC][k][j][i],d->Vc[TRC+1][k][j][i])*CONST_amu);
      
      q_m =  n_H*n_H*scrh/unit_q - LamHeat(i-2); //code units
     
     /* 
      if (g_stepNumber == 0){
      printf("%f  %f   \%f\n", n_H*n_H*scrh/unit_q,LamHeat(i-2),q_m);
     }
      */
      
      dq = d->Vc[PRS][k][j][i]*(1. - 1./(1. + q_m*dtsub*(g_gamma-1.)/d->Vc[PRS][k][j][i]));
      
      if (x1[i]> xbeg+noch && x1[i]< xend-noch){
        d->Vc[PRS][k][j][i] -= dq;
      }
      
      
      
    /*  if (cnt==0 && i==(int)NX1/2){
      printf("%0.100lf | %ld\n", d->Vc[PRS][k][j][i],g_stepNumber);}
     */ 
      
      
      //d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][i]/(1. + q_m*dtsub*(g_gamma-1.)/d->Vc[PRS][k][j][i]);
      //q_p += dq*GXYZ->dV[k][j][i]/(g_gamma-1.)/dtsub;
      //vol += GXYZ->dV[k][j][i];
    } /* -- end loop on points -- */
   /* 
    sendarray[0] = q_p; sendarray[1] = vol;
    #ifdef PARALLEL
    MPI_Allreduce (sendarray, recvarray, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #else
    recvarray[0] = sendarray[0]; recvarray[1] = sendarray[1];
    #endif
    q_p = recvarray[0]/recvarray[1]; //heating rate = vol. avg. cooling rate in code units

    DOM_LOOP(k,j,i){
      d->Vc[PRS][k][j][i] += (g_gamma-1.)*dtsub*q_p;
    }
  */
  } /* -- for loop cnt -- */
}
