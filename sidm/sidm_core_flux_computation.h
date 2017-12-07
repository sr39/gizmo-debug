/* here is where we call the core of the SIDM calculation for DM particle-particle interactions */
/* check if target particle is an SIDM candidate */
{
    if((1 << ptype) & (DM_SIDM))
    {
        /* ok, now check if neighbor particle is also SIDM-active */
        if((1 << P[no].Type) & (DM_SIDM))
        {
            /* ok, now check against self-interactions */
            if(targetID != P[no].ID)
            {
                sidm_tstart = my_second();
                r = sqrt(r2);
#if defined(ADAPTIVE_GRAVSOFT_FORALL)
                double h_si = DMAX(targeth_si, All.SIDMSmoothingFactor * DMAX(PPP[no].AGS_Hsml,All.ForceSoftening[P[no].Type]));
#else
                double h_si = DMAX(targeth_si, All.SIDMSmoothingFactor * All.ForceSoftening[P[no].Type]);
#endif
                if(r < 2.0*h_si)
                {
                    double prob = prob_of_interaction(P[no].Mass, r, h_si, targetVel, P[no].Vel, targetdt_step);
                    if(prob > max_prob) max_prob = prob;
                        
                        if(prob > 0.2)
                        {
                            if(targetdt_step_sidm == 0 ||
                               prob_of_interaction(P[no].Mass, r, h_si, targetVel, P[no].Vel, targetdt_step_sidm) > 0.2)
                            {
                                targetdt_step_sidm = targetdt_step;
                                double prob_tmp = prob;
                                while(prob_tmp > 0.2)
                                {
                                    targetdt_step_sidm /= 2;
                                    prob_tmp = prob_of_interaction(P[no].Mass, r, h_si, targetVel, P[no].Vel, targetdt_step_sidm);
                                }
                            }
                        } // if(prob > 0.2)
                    
                    if (gsl_rng_uniform(random_generator) < prob)
                    {
                        if(check_interaction_table(targetID, P[no].ID) == 0)
                        {
                            double kick_target[3], kick_no[3];
                            calculate_interact_kick(targetVel, P[no].Vel, kick_target, kick_no);
                            sidm_kick_x += kick_target[0];
                            sidm_kick_y += kick_target[1];
                            sidm_kick_z += kick_target[2];
                            int k;
                            for (k = 0; k < 3 ; k++)
                                P[no].Vel[k] += kick_no[k];
                                si_count++;
                            P[no].NInteractions++;
                            update_interaction_table(targetID, P[no].ID);
                        }  // if(check_interaction_table(targetID, P[no].ID) == 0)
                    } // if(prob for kick satisfied)
                } // if(r < 2.0*h_si)
                sidm_tend = my_second();
                sidm_tscatter += timediff(sidm_tstart, sidm_tend);
            } // if(targetID != P[no].ID)
        } // if((1 << P[no].Type) & (DM_SIDM))
    } // if((1 << ptype) & (DM_SIDM))
}

