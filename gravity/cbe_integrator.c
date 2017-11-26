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
void do_cbe_flux_computation(double *moments, double vface_dot_A, double *Area, double *fluxes);
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


/* this computes the actual single-sided fluxes at the face, integrating over a distribution function to use the moments */
void do_cbe_flux_computation(double *moments, double vface_dot_A, double *Area, double *fluxes)
{
    if(moments[0] <= 0) // no mass, no flux
    {
        int k; for(k=0;k<10;k++) {fluxes[k]=0;}
        return;
    }
    // couple dot-products must be pre-computed for fluxes //
    double m_inv = 1. / moments[0]; // need for weighting, below [e.g. v_x = moments[1] / moments[0]]
    double v[3]; v[0] = m_inv*moments[1]; v[1] = m_inv*moments[2]; v[2] = m_inv*moments[3]; // get velocities
    double v_dot_A = v[0]*Area[0] + v[1]*Area[1] + v[2]*Area[2]; // v_alpha . A_face
    double S[6]; // dispersion part of stress tensor (need to subtract mean-v parts
    S[0] = m_inv*moments[4] - v[0]*v[0]; // xx
    S[1] = m_inv*moments[5] - v[1]*v[1]; // yy
    S[2] = m_inv*moments[6] - v[2]*v[2]; // zz
    S[4] = m_inv*moments[7] - v[0]*v[1]; // xy
    S[5] = m_inv*moments[8] - v[0]*v[2]; // xz
    S[6] = m_inv*moments[9] - v[1]*v[2]; // yz
    double S_dot_A[3]; // what we actually use is this dotted into the face
    S_dot_A[0] = (S[0]*Area[0] + S[4]*Area[1] + S[5]*Area[2]) * moments[0]; // (S_alpha . A_face)_x * mass
    S_dot_A[1] = (S[4]*Area[0] + S[1]*Area[1] + S[6]*Area[2]) * moments[0]; // (S_alpha . A_face)_y * mass
    S_dot_A[2] = (S[5]*Area[0] + S[6]*Area[1] + S[2]*Area[2]) * moments[0]; // (S_alpha . A_face)_z * mass
    fluxes[0] = (v_dot_A - vface_dot_A) * moments[0]; // calculate and assign mass flux
    int k; for(k=1;k<10;k++) {fluxes[k] =  (m_inv * moments[k]) * fluxes[0];} // specific flux carried just by mass flux
    fluxes[1] += S_dot_A[0]; // add momentum flux from stress tensor - x
    fluxes[2] += S_dot_A[1]; // add momentum flux from stress tensor - y
    fluxes[3] += S_dot_A[2]; // add momentum flux from stress tensor - z
    fluxes[4] += 2. * v[0]*S_dot_A[0]; // add stress flux from stress tensor -- xx
    fluxes[5] += 2. * v[1]*S_dot_A[1]; // add stress flux from stress tensor -- yy
    fluxes[6] += 2. * v[2]*S_dot_A[2]; // add stress flux from stress tensor -- zz
    fluxes[7] += v[0]*S_dot_A[1] + v[1]*S_dot_A[0]; // add stress flux from stress tensor -- xy
    fluxes[8] += v[0]*S_dot_A[2] + v[2]*S_dot_A[0]; // add stress flux from stress tensor -- xz
    fluxes[9] += v[1]*S_dot_A[2] + v[2]*S_dot_A[1]; // add stress flux from stress tensor -- yz
    return;
}


#endif
