#include "pluto.h"
#include "omega.h"
#include "eigenmode.h"
#include "background.h"
#include "user_globals.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
   double kw, m, c, mue=1.17, mu=0.62;
   int iso;
   
   g_gamma = 5./3.; g_minCoolingTemp = 1.e6;
   kw = g_kw;  // Defined in user_globals.h
   iso = g_iso;
   
   m = g_m;
   c = g_cz;

/*
   us[VX1] = 0.0;
   us[RHO] = 0.1*(CONST_amu*mu/UNIT_DENSITY)*(1. + g_inputParam[amp]*sin(kw*x1));
   us[PRS] = 0.1*CONST_kB*0.78*1.16e7/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
   us[TRC] = 0.01*(sin(kw*x1)) + m*x1 + c; //this is metallicity in solar units; Z/Zsun and not Z!
*/

   double rho0 = rhobg(x1);
   double del_rho = 0.1*(CONST_amu*mu/UNIT_DENSITY)*g_inputParam[amp];
   us[RHO] = rho0 + del_rho*sin(kw*x1);
   //printf("rho done\n");
   
   double del_vr = delvr(kw,x1,del_rho);
   double del_vi = delvi(kw,x1,del_rho);
   double v0 = vbg(x1);
   
   us[VX1] = v0 + mag(del_vr,del_vi) * sin( kw*x1 + arg(del_vr,del_vi) ) ;
   //printf("p done\n");
   
   double del_pr = delpr(kw,x1,del_rho);
   double del_pi = delpi(kw,x1,del_rho);
   double p0 = pbg(x1);
   
   us[PRS] = p0 + mag(del_pr,del_pi) * sin( kw*x1 + arg(del_pr,del_pi) );
   //printf("v done\n");
   
   double del_zr = 0.0;
   double del_zi = delzi(kw,x1,del_rho);
   double z0 = zbg(x1);

   us[TRC] = z0 + mag(del_zr,del_zi) * sin( kw*x1 + arg(del_zr,del_zi) );    
   //printf("z done\n");
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{
  int k, j, i;
  double g_mass, g_TE, g_KE1, g_KE2, g_KE3, g_mom1, g_mom2, g_mom3, dvol;
  FILE *hist_file;

  #ifdef PARALLEL
  if (prank==0)
  #endif
  if (g_stepNumber==0) {
    hist_file = fopen ("pluto_hst.out", "w");
    fprintf(hist_file,"#[1]time [2]g_dt [3]mass [4]TE [5]KE1 [6]KE2 [7]KE3 [8]MOM1 [9]MOM2 [10]MOM3\n ");
  }
  else hist_file = fopen ("pluto_hst.out", "a");

  g_mass=0.0; g_TE=0.0; g_KE1=0.0; g_KE2=0.0; g_KE3=0.0;
  g_mom1=0.0; g_mom2=0.0; g_mom3=0.0;

  DOM_LOOP(k,j,i){
        dvol = grid->dV[k][j][i];
        g_mass += d->Vc[RHO][k][j][i]*dvol;
        g_TE += d->Vc[PRS][k][j][i]*dvol/(g_gamma-1.0);
        EXPAND( g_KE1 += 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*d->Vc[VX1][k][j][i]*dvol;  ,
                g_KE2 += 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX2][k][j][i]*d->Vc[VX2][k][j][i]*dvol;  ,
                g_KE3 += 0.5*d->Vc[RHO][k][j][i]*d->Vc[VX3][k][j][i]*d->Vc[VX3][k][j][i]*dvol; )
        EXPAND( g_mom1 += d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*dvol;  ,
                g_mom2 += d->Vc[RHO][k][j][i]*d->Vc[VX2][k][j][i]*dvol;  ,
                g_mom3 += d->Vc[RHO][k][j][i]*d->Vc[VX3][k][j][i]*dvol; )
  }

  #ifdef PARALLEL
   sendarray[0]=g_mass; sendarray[1]=g_TE; sendarray[2]=g_KE1; sendarray[3]=g_KE2;
   sendarray[4]=g_KE3; sendarray[5]=g_mom1; sendarray[6]=g_mom2; sendarray[7]=g_mom3;
   MPI_Reduce (sendarray, recvarray, 8, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
   if (prank == 0){
     g_mass=recvarray[0]; g_TE=recvarray[1]; g_KE1=recvarray[2]; g_KE2=recvarray[3];
     g_KE3=recvarray[4]; g_mom1=recvarray[5]; g_mom2=recvarray[6]; g_mom3=recvarray[7];
  #endif
    fprintf(hist_file,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n ", g_time, g_dt, g_mass, g_TE, g_KE1, g_KE2, g_KE3, g_mom1, g_mom2, g_mom3);

    fclose(hist_file);
  #ifdef PARALLEL
    }
  #endif  

}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   ib, ie, ieg, ibg, i, j, k, l;
  double T0;
  
  if (side ==0){
    TOT_LOOP(k,j,i){
        T0  = d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]*KELVIN*MMWt_mu(d->Vc[TRC][k][j][i]); //in cgs units
                
        if (T0 > g_ceil){
            d->Vc[PRS][k][j][i] = g_ceil*d->Vc[RHO][k][j][i]/(KELVIN*MMWt_mu(d->Vc[TRC][k][j][i]));
        }
        if (T0 < g_minCoolingTemp){
            d->Vc[PRS][k][j][i] = g_minCoolingTemp*d->Vc[RHO][k][j][i]/(KELVIN*MMWt_mu(d->Vc[TRC][k][j][i]));   
        }
      }
    }
  
  if (side == X1_BEG){
    BOX_LOOP(box,k,j,i){ 
      
      ib = IBEG;
      ie = IEND;

      ieg = ib + NX1;
      ibg = ie - NX1;     
      
      for(l=0;l<2;l++)
      {
        d->Vc[RHO][k][j][ibg-l] = d->Vc[RHO][k][j][ib+l];
        d->Vc[PRS][k][j][ibg-l] = d->Vc[PRS][k][j][ib+l];
        d->Vc[VX1][k][j][ibg-l] = d->Vc[VX1][k][j][ib+l];
        d->Vc[TRC][k][j][ibg-l] = d->Vc[TRC][k][j][ib+l];
        
	    d->Vc[RHO][k][j][ieg+l] = d->Vc[RHO][k][j][ie-l];
        d->Vc[PRS][k][j][ieg+l] = d->Vc[PRS][k][j][ie-l];
        d->Vc[VX1][k][j][ieg+l] = d->Vc[VX1][k][j][ie-l]; 
        d->Vc[TRC][k][j][ieg+l] = d->Vc[TRC][k][j][ie-l];

      }   

   }
  }
  /*else if (side == X1_END){
    BOX_LOOP(box,k,j,i){
     
      int ip = i - NX1;// - 1;
      d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][ip];
      d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][ip];
      d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][ip];
      d->Vc[TRC][k][j][i] = d->Vc[TRC][k][j][ip] - (m/x1[ip]) + (m/x1[i]);

      fi = fopen("outputend.txt","a");
      fprintf(fi,"%d %d\t",i,ip);

      int ib = i - 2;
      ip = ib - NX1;

      d->Vc[RHO][k][j][ib] = d->Vc[RHO][k][j][ip];
      d->Vc[PRS][k][j][ib] = d->Vc[PRS][k][j][ip];
      d->Vc[VX1][k][j][ib] = d->Vc[VX1][k][j][ip];
      d->Vc[TRC][k][j][ib] = d->Vc[TRC][k][j][ip] - (m/x1[ip]) + (m/x1[ib]);

      fprintf(fi,"%d %d\t",ib,ip);
      fclose(fi);

      //printf("end:%d  %f \nend:%d  %f \n\n",i,d->Vc[RHO][k][j][i],ip,d->Vc[RHO][k][j][ip]);
   }
  } */

}

