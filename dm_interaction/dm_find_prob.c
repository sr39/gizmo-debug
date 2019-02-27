//
//  dm_find_prob.c
//  ダークマターとガス粒子が衝突する確率を与える。
//
//  Created by 市橋 on 2016/11/30.
//
//

#include "dm_find_prob.h"
#include "dm_interaction.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"




#ifdef DM_BARYON_INTERACTION
void do_sph_interaction_kick(int i, integertime tstart, integertime tend, double dt_entr)
{
   int  k;
if(SphP[i].dm_DtInternalEnergy != 0){
printf("Sph[i].dtinternal = %e\n",SphP[i].dm_DtInternalEnergy);
printf("Sph[i].internal = %e\n",SphP[i].InternalEnergy);
printf("SphP[i].vel = %e\n",SphP[i].VelPred[1]);
                                    }
if(SphP[i].baryon_dtVel[1] != 0){
printf("Sph[i].dtvel[1] = %e\n",SphP[i].baryon_dtVel[1]);
                               }   
   SphP[i].InternalEnergy += SphP[i].dm_DtInternalEnergy * dt_entr;

   for(k = 0; k < 3; k++)
    {
       P[i].Vel[k] += SphP[i].baryon_dtVel[k] * dt_entr;
       
    }

}

void count_interaction(void)
{
 int i, count0 = 0, count1 = 0;
 for(i = 0; i < NumPart; i++)
  {
  if(P[i].Type == 0)
   {
    count0 += SphP[i].count0;
   }else{   
    count1 += P[i].count1;
   }
  }
   if(ThisTask==0)
   {
     printf("interaction_gas = %d\n", count0);
     printf("interaction_dm = %d\n", count1);
   }
}

void count_reset(void)
{
  int i;
  for(i = 0; i < NumPart; i++)
  {
   P[i].count1 = 0;
  }
}
#endif
