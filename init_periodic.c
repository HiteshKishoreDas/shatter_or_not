#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
   double kw, m, c, mue=1.17, mu=0.62;

   g_gamma = 5./3.; g_minCoolingTemp = 1.e6;
   kw = 2.*CONST_PI;
   m = 1.;
   c = 2.;

   us[VX1] = 0.0;
   //us[RHO] = 0.1*(CONST_amu*mue/UNIT_DENSITY)*(1. + g_inputParam[amp]*sin(w*,d->Vc[RHO][k][j][i],i+NX1+1,d->Vc[RHO][k][j][i+NX1+1]);x1));
   //us[PRS] = (0.1*mue/mu)*CONST_kB*0.78*1.16e7/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY); //2010 parameters
   us[RHO] = 0.1*(CONST_amu*mu/UNIT_DENSITY)*(1. + g_inputParam[amp]*sin(kw*x1));
   us[PRS] = 0.1*CONST_kB*0.78*1.16e7/(UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY);
   us[TRC] = 0.1*(sin(kw*x1)) ; //this is metallicity in solar units; Z/Zsun and not Z!
   //us[TRC] = 0.3*(1. + 0.*g_inputParam[amp]*sin(kw*x1)); //tracer stands for metallicity in this setup

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
  int   ib, ie, ipb, ipe, l, i, j, k, i1, j1;
  double  *x1,m,c;
  FILE *fi;
  //int cnt0=0,cnt1=0;
  m = 1.; 
  c = 2.;

  x1 = grid->x[IDIR];
  //printf("length: %f %f %f \n",x1[2049],x1[2050],x1[2051]);
  int count =0;
  if (side == X1_BEG){
    BOX_LOOP(box,k,j,i){ 
      
      ib = IBEG;
      ie = IEND;

      ipb = ib + NX1;
      ipe = ie - NX1;

     // printf("%f, %f | %f , %f || %f, %f | %f, %f\n", d->Vc[RHO][k][j][ib-2], d->Vc[RHO][k][j][ib-1], d->Vc[RHO][k][j][ib] , d->Vc[RHO][k][j][ib+1], d->Vc[RHO][k][j][ie-1], d->Vc[RHO][k][j][ie],d->Vc[RHO][k][j][ie+1],d->Vc[RHO][k][j][ib+2]);
      
      for(l=0;l<2;l++)
      {
        d->Vc[RHO][k][j][ipb+l] = d->Vc[RHO][k][j][ib+l];
        d->Vc[PRS][k][j][ipb+l] = d->Vc[PRS][k][j][ib+l];
        d->Vc[VX1][k][j][ipb+l] = d->Vc[VX1][k][j][ib+l];
        d->Vc[TRC][k][j][ipb+l] = d->Vc[TRC][k][j][ib+l];// + (c - m*x1[ipb+l]) - (c - m*x1[ib+l]);
        //printf("%d\n", l); 
	d->Vc[RHO][k][j][ipe-l] = d->Vc[RHO][k][j][ie-l];
        d->Vc[PRS][k][j][ipe-l] = d->Vc[PRS][k][j][ie-l];
        d->Vc[VX1][k][j][ipe-l] = d->Vc[VX1][k][j][ie-l];
        d->Vc[TRC][k][j][ipe-l] = d->Vc[TRC][k][j][ie-l];// + (c - m*x1[ipe-l]) - (c - m*x1[ie-l]);

       // printf("%f , %f\n",d->Vc[RHO][k][j][ib+l] , d->Vc[RHO][k][j][ipb+l]);
        //`printf("%f , %f\n",d->Vc[RHO][k][j][ie-l] , d->Vc[RHO][k][j][ipe-l]);
      }   
            
       
      double trb0 = d->Vc[RHO][k][j][ib-2];// - (c - m*x1[ib-2]);
      double trb1 = d->Vc[RHO][k][j][ib-1];// - (c - m*x1[ib-1]);
      double trb2 = d->Vc[RHO][k][j][ib];// - (c - m*x1[ib]);
      double trb3 = d->Vc[RHO][k][j][ib+1];// - (c - m*x1[ib+1]);

      double tre0 = d->Vc[RHO][k][j][ie-1];// - (c - m*x1[ie-1]);
      double tre1 = d->Vc[RHO][k][j][ie];// - (c - m*x1[ie]);
      double tre2 = d->Vc[RHO][k][j][ie+1];// - (c - m*x1[ie+1]);
      double tre3 = d->Vc[RHO][k][j][ie+2];// - (c - m*x1[ie+2]);
       
      printf("%f, %f | %f , %f || %f, %f | %f, %f\n", trb0 , trb1, trb2, trb3, tre0, tre1, tre2, tre3);

      fi = fopen("output.txt","a");
      fprintf(fi,"%f %f %f %f\n",trb0,trb1,tre2,tre3);

      fclose(fi);

      
      /*
      int ip = i + NX1;// +1;
      d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][ip];
      d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][ip];
      d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][ip];
      d->Vc[TRC][k][j][i] = d->Vc[TRC][k][j][ip] - (m/x1[ip]) + (m/x1[i]);
     
      
      
      //fi = fopen("outputbeg.txt","a");
      //fprintf(fi,"%d %d\t ",i,ip);
      
      int ib = IBEG;
      int ip = IEND;
      
      d->Vc[RHO][k][j][ib] = d->Vc[RHO][k][j][ip];
      d->Vc[PRS][k][j][ib] = d->Vc[PRS][k][j][ip];
      d->Vc[VX1][k][j][ib] = d->Vc[VX1][k][j][ip];
      d->Vc[TRC][k][j][ib] = d->Vc[TRC][k][j][ip] - (m/x1[ip]) + (m/x1[ib]);
      */
      /*int l;  
      
      if (cnt0==0 || cnt1==0){
      if (g_stepNumber == 0 || g_stepNumber == 500)
      {
        for (l=0;l<10;l++)
        fprintf(fi,"%ld\t%d\t%f\n",g_stepNumber,ib+l,d->Vc[RHO][k][j][ib+l]);
        if(g_stepNumber==0) cnt0=1;
        else cnt1 = 1;
      }}

      fclose(fi);
      */
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

