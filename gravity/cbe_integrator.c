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

#ifdef CBE_INTEGRATOR



// moment ordering convention: 0, x, y, z, xx, yy, zz, xy, xz, yz
//                             0, 1, 2, 3,  4,  5,  6,  7,  8,  9


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
// we evolve conserved quantities directly (in physical units): minor conversion in fluxes required later //
void do_cbe_drift_kick(int i, double dt)
{
    int j, k;
    double moment[10]={0}, dmoment[10]={0}, minv=1./P[i].Mass, v0[3]={0};
    // evaluate total fluxes //
    for(j=0;j<CBE_INTEGRATOR_NBASIS;j++)
    {
        for(k=0;k<10;k++)
        {
            moment[k] += P[i].CBE_basis_moments[j][k];
            dmoment[k] += dt*P[i].CBE_basis_moments_dt[j][k];
        }
    }
#ifdef CBE_DEBUG
    // total mass flux should vanish identically -- check this //
    if(ThisTask==0) {printf("Total mass flux == %g (Task=%d i=%d)\n",dmoment[0],ThisTask,i);}
    if(ThisTask==0) {printf("Momentum flux == %g/%g/%g (Task=%d i=%d)\n",dmoment[1],dmoment[2],dmoment[3],ThisTask,i);}
    if(fabs(minv*dmoment[0]) > 0.05) endrun(91918282);
#endif
    // define the current velocity, force-sync update to match it //
    for(k=0;k<3;k++) {v0[k] = P[i].Vel[k] / All.cf_atime;} // physical units //
    double biggest_dm = 1.e10;
    for(j=0;j<CBE_INTEGRATOR_NBASIS;j++)
    {
        double q = (dt*P[i].CBE_basis_moments_dt[j][0]) / (P[i].CBE_basis_moments[j][0] * (1.+minv*dmoment[0]));
        if(!isnan(q)) {if(q < biggest_dm) {biggest_dm=q;}}
    }
    double nfac = 1; // normalization factor for fluxes below //
    double threshold_dm;
    threshold_dm = -0.9; // maximum allowed fractional change in m //
#ifndef CBE_DEBUG
    if(biggest_dm < threshold_dm) {nfac = threshold_dm/biggest_dm;} // re-normalize flux so it doesn't overshoot //
#endif
    // ok now do the actual update //
    for(j=0;j<CBE_INTEGRATOR_NBASIS;j++)
    {
        P[i].CBE_basis_moments[j][0] += nfac * (dt*P[i].CBE_basis_moments_dt[j][0] - P[i].CBE_basis_moments[j][0]*minv*dmoment[0]); // update mass (strictly ensuring total mass matches updated particle)
        for(k=1;k<4;k++) {P[i].CBE_basis_moments[j][k] += nfac * (dt*P[i].CBE_basis_moments_dt[j][k] - P[i].CBE_basis_moments[j][0]*(minv*dmoment[k]-v0[k-1]));} // update momentum (strictly ensuring total momentum matches updated particle)
        for(k=4;k<10;k++) {P[i].CBE_basis_moments[j][k] += nfac * (dt*P[i].CBE_basis_moments_dt[j][k]);} // pure dispersion, no re-normalization here
#ifdef CBE_DEBUG
        // check against negatives //
        if(P[i].CBE_basis_moments[j][0]<0 || P[i].CBE_basis_moments[j][4]<0 || P[i].CBE_basis_moments[j][5]<0 || P[i].CBE_basis_moments[j][6]<0)
        {
            printf("ZZa: Task=%d i=%d j=%d   : m=%g v=%g/%g/%g Sxx,yy,zz=%g/%g/%g Sxy,xz,yz=%g/%g/%g \n",ThisTask,i,j,
                   P[i].CBE_basis_moments[j][0],P[i].CBE_basis_moments[j][1],P[i].CBE_basis_moments[j][2],P[i].CBE_basis_moments[j][3],
                   P[i].CBE_basis_moments[j][4],P[i].CBE_basis_moments[j][5],P[i].CBE_basis_moments[j][6],P[i].CBE_basis_moments[j][7],
                   P[i].CBE_basis_moments[j][8],P[i].CBE_basis_moments[j][9]);
            printf("ZZb: Task=%d i=%d j=%d dt: m=%g v=%g/%g/%g Sxx,yy,zz=%g/%g/%g Sxy,xz,yz=%g/%g/%g \n",ThisTask,i,j,
                   P[i].CBE_basis_moments_dt[j][0],P[i].CBE_basis_moments_dt[j][1],P[i].CBE_basis_moments_dt[j][2],P[i].CBE_basis_moments_dt[j][3],
                   P[i].CBE_basis_moments_dt[j][4],P[i].CBE_basis_moments_dt[j][5],P[i].CBE_basis_moments_dt[j][6],P[i].CBE_basis_moments_dt[j][7],
                   P[i].CBE_basis_moments_dt[j][8],P[i].CBE_basis_moments_dt[j][9]);
            fflush(stdout);
        }
#endif
    }
    return;
}




/* routine to invert the NV_T matrix after neighbor pass */
double do_cbe_nvt_inversion_for_faces(int i)
{
    MyFloat NV_T[3][3]; int j,k;
    for(j=0;j<3;j++) {for(k=0;k<3;k++) {NV_T[j][k]=P[i].NV_T[j][k];}} // initialize matrix to be inverted //
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
    double S[6]; // dispersion part of stress tensor (need to subtract mean-v parts if not doing so in pre-step)
    // note that we are actually evolving S, although we will compute the -flux- of T, the conserved quantity //
    S[0] = m_inv*moments[4];// - v[0]*v[0]; // xx
    S[1] = m_inv*moments[5];// - v[1]*v[1]; // yy
    S[2] = m_inv*moments[6];// - v[2]*v[2]; // zz
    S[3] = m_inv*moments[7];// - v[0]*v[1]; // xy
    S[4] = m_inv*moments[8];// - v[0]*v[2]; // xz
    S[5] = m_inv*moments[9];// - v[1]*v[2]; // yz
    double S_dot_A[3]; // what we actually use is this dotted into the face
    S_dot_A[0] = (S[0]*Area[0] + S[3]*Area[1] + S[4]*Area[2]) * moments[0]; // (S_alpha . A_face)_x * mass
    S_dot_A[1] = (S[3]*Area[0] + S[1]*Area[1] + S[5]*Area[2]) * moments[0]; // (S_alpha . A_face)_y * mass
    S_dot_A[2] = (S[4]*Area[0] + S[5]*Area[1] + S[2]*Area[2]) * moments[0]; // (S_alpha . A_face)_z * mass
    fluxes[0] = (v_dot_A - vface_dot_A) * moments[0]; // calculate and assign mass flux
    int k; for(k=1;k<10;k++) {fluxes[k] =  (m_inv * moments[k]) * fluxes[0];} // specific flux carried just by mass flux
    fluxes[1] += S_dot_A[0]; // add momentum flux from stress tensor - x
    fluxes[2] += S_dot_A[1]; // add momentum flux from stress tensor - y
    fluxes[3] += S_dot_A[2]; // add momentum flux from stress tensor - z
    fluxes[4] += 2.*v[0]*S_dot_A[0] + fluxes[0]*v[0]*v[0]; // add stress flux from stress tensor -- xx
    fluxes[5] += 2.*v[1]*S_dot_A[1] + fluxes[0]*v[1]*v[1]; // add stress flux from stress tensor -- yy
    fluxes[6] += 2.*v[2]*S_dot_A[2] + fluxes[0]*v[2]*v[2]; // add stress flux from stress tensor -- zz
    fluxes[7] += v[0]*S_dot_A[1] + v[1]*S_dot_A[0] + fluxes[0]*v[0]*v[1]; // add stress flux from stress tensor -- xy
    fluxes[8] += v[0]*S_dot_A[2] + v[2]*S_dot_A[0] + fluxes[0]*v[0]*v[2]; // add stress flux from stress tensor -- xz
    fluxes[9] += v[1]*S_dot_A[2] + v[2]*S_dot_A[1] + fluxes[0]*v[1]*v[2]; // add stress flux from stress tensor -- yz
    return;
}





/* this routine contains operations which are needed after the main forcetree-walk loop (where the CBE integration terms are calculated).
     we shift the net momentum flux into the GravAccel vector so that the tree and everything else behaves correctly, velocities
     drift, etc, all as they should. The 'residual' terms are then saved, which can kicked separately from the main particle kick.
*/
void do_postgravity_cbe_calcs(int i)
{
    int j,k;
    double (*mom)[10] = P[i].CBE_basis_moments;
    double (*dmom)[10] = P[i].CBE_basis_moments_dt;
    double dmom_tot[10]={0}, m_inv = 1./P[i].Mass;
    for(j=0;j<CBE_INTEGRATOR_NBASIS;j++) {for(k=0;k<10;k++) {dmom_tot[k] += dmom[j][k];}} // total change for each moment
#ifdef CBE_DEBUG
    if(ThisTask==0) {printf("mom=%g/%g/%g/%g \n",mom[0][3],mom[2][7],mom[1][5],mom[0][8]);}
    /* total mass change should be zero to floating-point accuracy, so don't need to worry about it, but check! */
    if(ThisTask==0) {printf("PG: Total mass flux == %g (Task=%d i=%d)\n",dmom_tot[0],ThisTask,i);}
    if(ThisTask==0) {printf("PG: Momentum flux == %g/%g/%g (Task=%d i=%d)\n",dmom_tot[1],dmom_tot[2],dmom_tot[3],ThisTask,i);}
#endif
    /* total momentum flux will be transferred */
    double dv0[3];
    for(k=0;k<3;k++)
    {
        dv0[k] = m_inv * dmom_tot[k]; // total acceleration
        P[i].GravAccel[k] += dv0[k] / All.cf_a2inv; // write as gravitational acceleration, convert to cosmological units
    }
    // now need to add that shift back into the momentum-change terms //
    for(j=0;j<CBE_INTEGRATOR_NBASIS;j++)
    {
        // first, shift the second-moment derivatives so they are dS, not dT, where T = S + v.v (outer product v.v),
        //   so dS = dT - (dv.v + v.dv) = dT - (dv.v + transpose[dv.v])
        double dS[6]={0};
        dS[0] = dmom[j][4] - m_inv * (dmom[j][1]*mom[j][1] + mom[j][1]*dmom[j][1]) + m_inv*m_inv * dmom[j][0] * mom[j][1]*mom[j][1]; // xx
        dS[1] = dmom[j][5] - m_inv * (dmom[j][2]*mom[j][2] + mom[j][2]*dmom[j][2]) + m_inv*m_inv * dmom[j][0] * mom[j][2]*mom[j][2]; // yy
        dS[2] = dmom[j][6] - m_inv * (dmom[j][3]*mom[j][3] + mom[j][3]*dmom[j][3]) + m_inv*m_inv * dmom[j][0] * mom[j][3]*mom[j][3]; // zz
        dS[3] = dmom[j][7] - m_inv * (dmom[j][1]*mom[j][2] + mom[j][2]*dmom[j][1]) + m_inv*m_inv * dmom[j][0] * mom[j][1]*mom[j][2]; // xy
        dS[4] = dmom[j][8] - m_inv * (dmom[j][1]*mom[j][3] + mom[j][3]*dmom[j][1]) + m_inv*m_inv * dmom[j][0] * mom[j][1]*mom[j][3]; // xz
        dS[5] = dmom[j][9] - m_inv * (dmom[j][2]*mom[j][3] + mom[j][4]*dmom[j][2]) + m_inv*m_inv * dmom[j][0] * mom[j][2]*mom[j][3]; // yz
        dmom[j][0] -= mom[j][0]*(m_inv * dmom_tot[0]); // re-ensure that this is zero to floating-point precision //
        for(k=1;k<4;k++) {dmom[j][k] -= mom[j][0]*dv0[k];} // shift the momentum flux, so now zero net 'residual' momentum flux //
        for(k=4;k<10;k++) {dmom[j][k] = dS[k-4];} // set the flux of the stress terms to the change in the dispersion alone //
    } // for(j=0;j<CBE_INTEGRATOR_NBASIS;j++)
    return;
}


#endif
