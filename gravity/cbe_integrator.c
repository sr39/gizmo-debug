#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"


/*! \file cbe_integrator.c
 *  \brief routines needed for CBE integrator implementation
 *         This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#define CBE_DEBUG

#ifdef CBE_INTEGRATOR
void do_cbe_initialization(void);
void do_cbe_drift_kick(int i, double dt);
double do_cbe_nvt_inversion_for_faces(int i);
#endif

// moment ordering convention: 0, x, y, z, xx, yy, zz, xy, xz, yz


#ifdef CBE_INTEGRATOR


/* variable initialization */
void do_cbe_initialization(void)
{
    int i,j,k;
    for(i=0;i<NumPart;i++)
    {
        for(j=0;j<CBE_INTEGRATOR_NBASIS;j++) {for(k=0;k<10;k++) {P[i].CBE_basis_moments_dt[j][k]=0;}} // no time derivatives //
        double mom_tot[10]={0};
        for(j=0;j<CBE_INTEGRATOR_NBASIS;j++)
        {
            for(k=0;k<10;k++)
            {
                // zeros will be problematic, instead initialize a random distribution //
                P[i].CBE_basis_moments[j][k] = 1.e-10 * P[i].Mass * get_random_number(P[i].ID + i + 343*ThisTask + 912*k + 781*j);
                if(k>0 && k<4) {P[i].CBE_basis_moments[j][k] *= P[i].Vel[k];}
                if(k>=4) {P[i].CBE_basis_moments[j][k] *= P[i].Vel[k]*P[i].Vel[k];}
#if (NUMDIMS==1)
                if((k!=0)&&(k!=1)&&(k!=4)) {P[i].CBE_basis_moments[j][k] = 0;}
#endif
#if (NUMDIMS==2)
                if((k==3)||(k==6)||(k==8)||(k==9)) {P[i].CBE_basis_moments[j][k] = 0;}
#endif
                mom_tot[k] += P[i].CBE_basis_moments[j][k];
            }
        }
        P[i].CBE_basis_moments[j][0] *= P[i].Mass / mom_tot[0]; // must normalize to correct mass
        for(k=0;k<3;k++) {P[i].CBE_basis_moments[j][k+1] -= mom_tot[k+1];} // ensure IC has zero momentum relative to bulk
    }
    return;
}



/* drift-kick updates to distribution functions */
// we evolve conserved quantities directly: minor conversion in fluxes required later //
void do_cbe_drift_kick(int i, double dt)
{
    int j, k;
    double dmoment[10]={0}, minv=1./P[i].Mass;
    // total mass flux should vanish identically //
    for(j=0;j<CBE_INTEGRATOR_NBASIS;j++) {dm_tot += dt * CBE_basis_moments_dt[j][0];} // mass flux
#ifdef CBE_DEBUG
    if(ThisTask==0) {printf("Total mass flux == %g (Task=%d i=%d)\n",dm_tot,ThisTask,i);}
#endif
    
    // loop over basis functions //
    for(j=0;j<CBE_INTEGRATOR_NBASIS;j++)
    {
        CBE_basis_moments[j][0] += dt*CBE_basis_moments_dt[j][0]; // update mass
        dmoment[0] += CBE_basis_moments[j][0];
        for(k=1;k<10;k++)
        {
            CBE_basis_moments[j][k] += dt*CBE_basis_moments_dt[j][k]); // momentum or energy update
            dmoment[k] += CBE_basis_moments[j][k] * minv; // sum total linear momentum (or energy), divided by total mass
        }
    }
#ifdef CBE_DEBUG
    if(fabs(dmoment[0]-P[i].Mass) > 0.05*P[i].Mass) endrun(91918282);
#endif
    for(k=0;k<3;k++) // ??? need cosmological units on vel below; does momentum re-zero shift energies? (it should) //
    {
        P[i].Vel[k] += dmoment[k+1]; // shift to new center-of-momentum
        if(P[i].Type==0) {SphP[i].VelPred[k] += dmoment[k+1]} // gas needs predictor-drift as well
        for(j=0;j<CBE_INTEGRATOR_NBASIS;j++) {CBE_basis_moments[j][k+1] -= P[i].Mass * dmoment[k+1];} // subtract back off, so now momentum should normalize to zero
    }
    return;
}




/* routine to invert the NV_T matrix after neighbor pass */
double do_cbe_nvt_inversion_for_faces(int i)
{
    MyDouble *NV_T = P[i].NV_T; // matrix to be inverted //
    int j,k;
    double Tinv[3][3], FrobNorm=0, FrobNorm_inv=0, detT=0;
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {Tinv[j][k]=0;}}
    /* fill in the missing elements of NV_T (it's symmetric, so we saved time not computing these directly) */
    NV_T[1][0]=NV_T[0][1]; NV_T[2][0]=NV_T[0][2]; NV_T[2][1]=NV_T[1][2];
    /* Also, we want to be able to calculate the condition number of the matrix to be inverted, since
     this will tell us how robust our procedure is (and let us know if we need to expand the neighbor number */
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {FrobNorm += NV_T[j][k]*NV_T[j][k];}}
#if (NUMDIMS==1) // 1-D case //
    detT = NV_T[0][0];
    if(detT!=0 && !isnan(detT)) {Tinv[0][0] = 1/detT}; /* only one non-trivial element in 1D! */
#endif
#if (NUMDIMS==2) // 2-D case //
    detT = NV_T[0][0]*NV_T[1][1] - NV_T[0][1]*NV_T[1][0];
    if((detT != 0)&&(!isnan(detT)))
    {
        Tinv[0][0] = NV_T[1][1] / detT; Tinv[0][1] = -NV_T[0][1] / detT;
        Tinv[1][0] = -NV_T[1][0] / detT; Tinv[1][1] = NV_T[0][0] / detT;
    }
#endif
#if (NUMDIMS==3) // 3-D case //
    detT = NV_T[0][0] * NV_T[1][1] * NV_T[2][2] + NV_T[0][1] * NV_T[1][2] * NV_T[2][0] +
    NV_T[0][2] * NV_T[1][0] * NV_T[2][1] - NV_T[0][2] * NV_T[1][1] * NV_T[2][0] -
    NV_T[0][1] * NV_T[1][0] * NV_T[2][2] - NV_T[0][0] * NV_T[1][2] * NV_T[2][1];
    /* check for zero determinant */
    if((detT != 0) && !isnan(detT))
    {
        Tinv[0][0] = (NV_T[1][1] * NV_T[2][2] - NV_T[1][2] * NV_T[2][1]) / detT;
        Tinv[0][1] = (NV_T[0][2] * NV_T[2][1] - NV_T[0][1] * NV_T[2][2]) / detT;
        Tinv[0][2] = (NV_T[0][1] * NV_T[1][2] - NV_T[0][2] * NV_T[1][1]) / detT;
        Tinv[1][0] = (NV_T[1][2] * NV_T[2][0] - NV_T[1][0] * NV_T[2][2]) / detT;
        Tinv[1][1] = (NV_T[0][0] * NV_T[2][2] - NV_T[0][2] * NV_T[2][0]) / detT;
        Tinv[1][2] = (NV_T[0][2] * NV_T[1][0] - NV_T[0][0] * NV_T[1][2]) / detT;
        Tinv[2][0] = (NV_T[1][0] * NV_T[2][1] - NV_T[1][1] * NV_T[2][0]) / detT;
        Tinv[2][1] = (NV_T[0][1] * NV_T[2][0] - NV_T[0][0] * NV_T[2][1]) / detT;
        Tinv[2][2] = (NV_T[0][0] * NV_T[1][1] - NV_T[0][1] * NV_T[1][0]) / detT;
    }
#endif
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {FrobNorm_inv += Tinv[j][k]*Tinv[j][k];}}
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {NV_T[j][k]=Tinv[j][k];}} // now NV_T holds the inverted matrix elements //
    double ConditionNumber = DMAX(sqrt(FrobNorm * FrobNorm_inv) / NUMDIMS, 1); // = sqrt( ||NV_T^-1||*||NV_T|| ) :: should be ~1 for a well-conditioned matrix //
#ifdef CBE_DEBUG
    if(ThisTask==0) {printf("Condition number == %g (Task=%d i=%d)\n",ConditionNumber,ThisTask,i);}
#endif
    return ConditionNumber;
}




#endif
