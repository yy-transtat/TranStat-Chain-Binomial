
double community_log_likelihood(COMMUNITY *community, double *par_effective)
{
  int i, j, h, k, l, m, n;
  int representative, converge, error=0;
  int start_day, stop_day, day_last_obs, min_size;
  int inf1, inf2;
  int r, rr, t, loop, id, index, skip, verify;
  int found, b_mode, p_mode, q_mode, u_mode;
  int n_par, n_b_mode, n_p_mode, n_q_mode, n_u_mode;
  int n_c2p_covariate, n_p2p_covariate, n_imm_covariate, n_pat_covariate;
  int n_covariate, n_time_ind_covariate, n_time_dep_covariate;
  int n_sus_p2p_covariate, n_inf_p2p_covariate, n_int_p2p_covariate; 
  int n_par_equiclass;

  double factor1, factor2, factor3;
  double e, f, ff, s, p_inf, Q, U;
  double covariate_effect, cum_e, day_L, L, LQ, log_L, my_log_L, my_L, log_L_all;
  
  double *par, sdf, pdf, pr;
  double logit_f, sum_prob;
  double temp, asym_effect;

  RISK *ptr_risk;
  PEOPLE *person, *member;
  RISK_CLASS *ptr_class;
  INTEGER_CHAIN *ptr_integer;


  /***********************************************************************************************************
   * n_b_mode: number of types/modes of community-to-person contact. 
   * n_p_mode: number of types/modes of person-to-person contact.
   * n_time_ind_covariate: number of time-independent covariates
   * n_time_dep_covariate: number of time-dependent covariates
   * n_covariate: total number of covariates
   * n_c2p_covariate: number of covariates of the susceptible person that can modify community-to-person risk
   * n_sus_p2p_covariate: number of covariates of the susceptible person that can modify person-to-person risk
   * n_inf_p2p_covariate: number of covariates of the infective person that can modify person-to-person risk
   * n_inf_p2p_covariate: number of interactions between covariates of the susceptible person
   *                      and the infective person that can modify person-to-person risk
   * n_p2p_covariate: total number of covariates that can modify person-to-person risk, which is the sum
   *                  of above three.
   * The user need to supply  n_time_ind_covariate, n_time_dep_covariate, n_c2p_covariate,
   * n_sus_p2p_covariate, n_inf_p2p_covariate and n_int_p2p_covariate.
   ***********************************************************************************************************/
  n_b_mode = cfg_pars.n_b_mode;
  n_p_mode = cfg_pars.n_p_mode;
  n_u_mode = cfg_pars.n_u_mode;
  n_q_mode = cfg_pars.n_q_mode;
  n_time_ind_covariate = cfg_pars.n_time_ind_covariate;
  n_time_dep_covariate = cfg_pars.n_time_dep_covariate;
  n_covariate = cfg_pars.n_covariate;
  n_c2p_covariate = cfg_pars.n_c2p_covariate;
  n_sus_p2p_covariate = cfg_pars.n_sus_p2p_covariate;
  n_inf_p2p_covariate = cfg_pars.n_inf_p2p_covariate;
  n_int_p2p_covariate = cfg_pars.n_int_p2p_covariate;
  n_p2p_covariate = cfg_pars.n_p2p_covariate;
  n_pat_covariate = cfg_pars.n_pat_covariate;
  n_imm_covariate = cfg_pars.n_imm_covariate;
  n_par = cfg_pars.n_par;
  n_par_equiclass = cfg_pars.n_par_equiclass;
  asym_effect = cfg_pars.asym_effect_est;
 

  make_1d_array_double(&par, n_par, 0.0);


  /************************************************************************************************************
   * f is person-to-person escape probability. e is the daily escape probability for a susceptible.
   * log_f_lb is the 1st derivative of log(f) with respect to lb.
   * log_f_lb_lb is the second derivative of log(f) with respect to lb.
   * Other derivative terms are similarly defined.
   *
   * e is the daily escape probability for a susceptible from all contacts.
   * cum_e and cum_log_e_* are cumulatives of e and log_e_* over a short period, 
   * from day_ill-max_latent to t, where day_ill-max_latent<=t<=day_ill-min_latent.
   * 
   * ee[t] is the escape probability for day t. 
   * cum_log_ee[t] = log(ee[1] * ee[2] * ... * ee[t])
   * 
   * temp_* variables are just auxiliary variables to help get cum_log_ee_* variables.
   * 
   * day_L is the likelihood for each infected person for the period 
   * from day_ill-max_latent+1 to t, where day_ill-max_latent<=t<=day_ill-min_latent.
   * log_day_L_* are derivatives of log(day_L). 
   *
   * L is the overall likelihood, and L_* and log_L_* are derivatives of L and log(L) respectively.
   * *********************************************************************************************************/
  //show_par_effective(par_effective, NULL);
  //make sure input parameter values are not out of scope
  for(i=0; i<n_par_equiclass; i++)
  {
     if(par_effective[i] < cfg_pars.lower_search_bound[i])
     {
        par_effective[i] = cfg_pars.lower_search_bound[i];
     }
     if(par_effective[i] > cfg_pars.upper_search_bound[i])
     {
        par_effective[i] = cfg_pars.upper_search_bound[i];
     }
     m = cfg_pars.par_equiclass[i].member[0] - 1;
     if(m < n_b_mode + n_p_mode)
     {
        if(par_effective[i] < logit(close_to_0))  par_effective[i] = logit(close_to_0);
        if(par_effective[i] > logit(close_to_05))  par_effective[i] = logit(close_to_05);
     }
     else if(m < n_b_mode + n_p_mode + n_u_mode + n_q_mode)
     {
        if(par_effective[i] < logit(close_to_0))  par_effective[i] = logit(close_to_0);
        if(par_effective[i] > logit(close_to_1))  par_effective[i] = logit(close_to_1);
     }
  }
  //show_par_effective(par_effective, NULL);
  for(i=0; i<n_par_equiclass; i++)
  {
     for(j=0; j<cfg_pars.par_equiclass[i].size; j++)
     {
        m = cfg_pars.par_equiclass[i].member[j] - 1;
        par[m] = par_effective[i];
     }
  }

  //if fixed parameters exist, need to set fixed parameters to prespecified values
  if(cfg_pars.n_par_fixed > 0)
  {
     for(i=0; i<cfg_pars.n_par_fixed; i++)
     {
        j = cfg_pars.par_fixed_id[i] - 1;
        par[j] =  cfg_pars.par_fixed_value[i];
     }
  }
  
  for(i=0; i<n_b_mode; i++)  {lb[i] = par[i]; b[i] = inv_logit(lb[i]);}  
  for(i=0; i<n_p_mode; i++)  {lp[i] = par[n_b_mode + i]; p[i] = inv_logit(lp[i]);}  
  for(i=0; i<n_u_mode; i++)  {lu[i] = par[n_b_mode + n_p_mode + i]; u[i] = inv_logit(lu[i]);}  
  for(i=0; i<n_q_mode; i++)  {lq[i] = par[n_b_mode + n_p_mode + n_u_mode + i]; q[i] = inv_logit(lq[i]);}  
  for(i=0; i<n_c2p_covariate; i++)  coeff_c2p[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode + i];  
  for(i=0; i<n_p2p_covariate; i++)  coeff_p2p[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode 
                                                 + n_c2p_covariate + i]; 
  for(i=0; i<n_pat_covariate; i++)  coeff_pat[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode 
                                                 + n_c2p_covariate + n_p2p_covariate + i]; 
  for(i=0; i<n_imm_covariate; i++)  coeff_imm[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode 
                                                 + n_c2p_covariate + n_p2p_covariate + n_pat_covariate + i]; 
  log_L_all = 0.0;
  if(cfg_pars.adjust_for_left_truncation == 1)  min_size = community->size - community->size_idx;
  else  min_size = community->size;
  //printf("min_size=%d\n", min_size);
  if(min_size > 0 && community->ignore == 0)
  {
     // index cases are now included in the risk classes.
     ptr_class = community->risk_class;
     while(ptr_class != NULL)
     {
        if(cfg_pars.common_contact_history_within_community == 1)
        {
           temp = 0.0;
           if(cfg_pars.adjust_for_left_truncation == 1 && community->earliest_idx_day_ill != MISSING)
              start_day = max(community->day_epi_start, community->earliest_idx_day_ill - cfg_pars.max_incubation + 1);
           else  start_day = community->day_epi_start;
           stop_day = community->day_epi_stop;
           //printf("start_day=%d  stop_day=%d  day_epi_stop=%d\n", start_day, stop_day, community->day_epi_stop);
           if(stop_day < start_day)
           {
              printf("community_log_likelihood, common contact history (line 164): community %d, stop day %d < start day %d\n", community->id, stop_day, start_day);
              exit(0);
           }
           
           for(t=start_day; t<=stop_day; t++)
           {
              r = t - start_day;
              rr = t - community->day_epi_start; //the time reference for risk history is community->day_epi_start, not start_day.
              e = 1.0;
              //printf("person %d: start_day=%d  t=%d  r=%d  rr=%d\n", person->id, start_day, t, r, rr);
              //calculate log_f_* variables related to the probability of escaping risk from common source.
              //update log_e_* variables accordingly
              
              ptr_risk = ptr_class->risk_history[rr].c2p_risk;
              while(ptr_risk != NULL)
              {
                 covariate_effect = 0.0;
                 for(k=0; k<n_c2p_covariate; k++)
                    covariate_effect += ptr_risk->covariate[k] * coeff_c2p[k];
                 covariate_effect += ptr_risk->offset;
                 
                 b_mode = ptr_risk->contact_mode;
                 logit_f = lb[b_mode] + covariate_effect;
                 
                 ff = inv_logit(logit_f);
                 f = 1.0 - ff;
                 
                 e *= ipow(f, ptr_risk->size);
                 ptr_risk = ptr_risk->next;
              }

              //calculate log_f_* variables related to the probability of escaping risk from infective people.
              //update log_e_* variables accordingly
              //The difference in infectiousness level between symptomatic and asymptomatic cases
              //is adjusted by one additional covariate that denote the symptom status of the infective person.
              ptr_risk = ptr_class->risk_history[rr].p2p_risk;
              while(ptr_risk != NULL)
              {
                 covariate_effect = 0.0;
                 for(k=0; k<n_p2p_covariate; k++)
                    covariate_effect += ptr_risk->covariate[k] * coeff_p2p[k];
                 covariate_effect += ptr_risk->offset;
                 
                 p_mode = ptr_risk->contact_mode;
                 logit_f = lp[p_mode] + covariate_effect;
                 s = ptr_risk->infective_prob;
                 //printf("id=%d  t=%d  s=%e  size=%d  p0=%e  p1=%e  ", 
                 //        person->id, t, s, ptr_risk->size, inv_logit(lp[p_mode]), inv_logit(logit_f));
                 //for(k=0; k<n_p2p_covariate; k++)  printf("x=%e  beta=%e  ", ptr_risk->covariate[k], coeff_p2p[k]);
                 //printf("\n");
                 if( s > 0)
                 {
                    ff = inv_logit(logit_f) * s * bipow(asym_effect, 1 - ptr_risk->symptom);
                    f = 1 - ff;
                    factor1 = ff / s;
                    
                    if(f <= 0.0 || f > 1.0)
                    {
                       printf("Full model: f=%e\n", f);
                       error = 1;
                       goto end;
                    }
                    e *= ipow(f, ptr_risk->size);
                 }
                 ptr_risk = ptr_risk->next;
              }

              //update log_ee_* variables, which are related the cumulative escape probabilities
              ee[r] = e;
              temp += log(ee[r]);
              cum_log_ee[r] = temp;
              //printf("%d: %e\n", t, ee[r]);
              //if(t == 310)  exit(0);
           } /* end of day t */
        }
        
        ptr_integer = ptr_class->member;
        while(ptr_integer != NULL)
        {
           person = people + ptr_integer->id;
           if(cfg_pars.adjust_for_left_truncation == 1 && person->idx == 1)  
           {
              ptr_integer = ptr_integer->next;
              continue;
           }   
           if(n_u_mode > 0)  u_mode = person->u_mode;
           if(n_q_mode > 0)  q_mode = person->q_mode;

           my_log_L = 0.0;
           // consider prior immunity levels
           if(n_q_mode > 0)
           {
              covariate_effect = 0.0;
              for(k=0; k<n_imm_covariate; k++)
                 covariate_effect += person->imm_covariate[k] * coeff_imm[k];
              Q = inv_logit(lq[q_mode] + covariate_effect);

              if( person->pre_immune == 1)
                    my_log_L += log(Q);
              else
                 my_log_L += log(1 - Q);
           }

           //show_p2p_risk(person, 1, 5);

           // if people in the comunity do not share the same exposure/risk history,
           // calculate individual-level exposure/risk history.
           if(cfg_pars.common_contact_history_within_community == 0 && person->pre_immune == 0)
           {
              temp = 0.0;
              if(cfg_pars.adjust_for_left_truncation == 1 && community->earliest_idx_day_ill != MISSING)
                 start_day = max(community->day_epi_start, community->earliest_idx_day_ill - cfg_pars.max_incubation + 1);
              else  start_day = community->day_epi_start;
              day_last_obs = min(community->day_epi_stop, person->day_exit);
              day_last_obs = min(community->day_last_followup, day_last_obs);
              stop_day = (person->infection == 1)? person->day_infection_upper:day_last_obs;
              if(stop_day < start_day)
              {
                 printf("community_log_likelihood (line 282): community %d, person %d, stop day %d < start day %d\n", community->id, person->id, stop_day, start_day);
                 exit(0);
              }
              for(t=start_day; t<=stop_day; t++)
              {
                 r = t - start_day;
                 rr = t - community->day_epi_start; //the time reference for risk history is community->day_epi_start, not start_day.

                 e = 1.0;
                 //calculate log_f_* variables related to the probability of escaping risk from community per day.
                 //update log_e_* variables accordingly
                 ptr_risk = person->risk_history[rr].c2p_risk;
                 while(ptr_risk != NULL)
                 {
                    covariate_effect = 0.0;
                    for(k=0; k<n_c2p_covariate; k++)
                       covariate_effect += ptr_risk->covariate[k] * coeff_c2p[k];
                    covariate_effect += ptr_risk->offset;
                    
                    b_mode = ptr_risk->contact_mode;
                    logit_f = lb[b_mode] + covariate_effect;

                    ff = inv_logit(logit_f);
                    f = 1.0 - ff;
                    e *= ipow(f, ptr_risk->size);
                    ptr_risk = ptr_risk->next;
                 }
                 //calculate log_f_* variables related to the probability of escaping risk from infective people.
                 //update log_e_* variables accordingly.
                 //The difference in infectiousness level between symptomatic and asymptomatic cases
                 //is adjusted by one additional covariate that denote the symptom status of the infective person.

                 ptr_risk = person->risk_history[rr].p2p_risk;
                 while(ptr_risk != NULL)
                 {
                    covariate_effect = 0.0;
                    for(k=0; k<n_p2p_covariate; k++)
                       covariate_effect += ptr_risk->covariate[k] * coeff_p2p[k];
                    covariate_effect += ptr_risk->offset;
                                
                    p_mode = ptr_risk->contact_mode;
                    logit_f = lp[p_mode] + covariate_effect;
                    s = ptr_risk->infective_prob;
                    if( s > 0)
                    {
                       ff = inv_logit(logit_f) * s * bipow(asym_effect, 1 - ptr_risk->symptom);
                       f = 1 - ff;

                       if(f <= 0.0 || f > 1.0)
                       {
                          //printf("Full model: f=%e\n", f);
                          error = 1;
                          goto end;
                       }
                       e *= ipow(f, ptr_risk->size);
                    }      
                    ptr_risk = ptr_risk->next;
                 }    /* if(person->p2p_contact_history[r].size > 0) */
                 ee[r] = e;

                 temp += log(ee[r]);
                 cum_log_ee[r] = temp;
              } /* end of day t */
           }
           if( person->pre_immune == 0)
           {
              if(cfg_pars.adjust_for_left_truncation == 1 && community->earliest_idx_day_ill != MISSING)
                    start_day = max(community->day_epi_start, community->earliest_idx_day_ill - cfg_pars.max_incubation + 1);
              else  start_day = community->day_epi_start;
              day_last_obs = min(community->day_epi_stop, person->day_exit);
              day_last_obs = min(community->day_last_followup, day_last_obs);
              stop_day = (person->infection == 1)? person->day_infection_upper:day_last_obs;
              if(stop_day < start_day)
              {
                 printf("community_log_likelihood (line 356): community %d, person %d, stop day %d < start day %d\n", community->id, person->id, stop_day, start_day);
                 printf("day_epi_stop=%d  day_exit=%d  day_last_followup=%d  day_infection_upper=%d\n", 
                         community->day_epi_stop, person->day_exit, community->day_last_followup, person->day_infection_upper);
                 exit(0);
              }
              else
              {
                 L = log_L = 0.0;

                 /******************************************************************************************
                  calculate likelihood for a infected case with symptoms 
                  ******************************************************************************************/
                    
                 if(person->infection == 1) /* if he was ill */
                 {
                    inf1 = person->day_infection_lower;
                    inf2 = person->day_infection_upper;
                    if(inf2 >= start_day)
                    {
                       // it could happen that inf1 < start_day < inf2
                       inf1 = max(inf1, start_day);
                       cum_e = 1.0;
                       for(t=inf1; t<=inf2; t++)
                       {
                          pdf = pdf_incubation[person->day_ill - cfg_pars.min_incubation - t];
                          r = t - start_day;
                          if(ee[r] < 0.0 || ee[r] > 1.0)
                          {
                             printf("community_loglikelihood (line 384), infected person id=%d: t=%d ee[r]=%e\n", person->id, t, ee[r]);
                             error = 1;
                             exit(0);
                          }
                             
                          if(ee[r] >= 0.0 && ee[r] < 1.0)
                          {
                             // if ee[r] is very close to 1, 1-ee[r] may just be set to 0 by the system.
                             // We found that -log(ee[r]) can solve this problem
                             p_inf = (fabs(log(ee[r]))<1e-6)? (-log(ee[r])) : (1.0 - ee[r]);
                             day_L = cum_e * p_inf;
                             if(n_u_mode > 0)
                             {
                                rr = t - community->day_epi_start; // track covariate for pathogenicity, for which timeline starts from community->day_epi_start.
                                                                   // Covariates affecting pathogenicity are from the day of infection 
                                covariate_effect = 0.0;
                                for(k=0; k<n_pat_covariate; k++)
                                   covariate_effect += person->pat_covariate[k][rr] * coeff_pat[k];
                                U = inv_logit(lu[u_mode] + covariate_effect);

                                pr = pdf * (person->symptom * U + (1 - person->symptom) * (1 - U));
                             }
                             else  pr = pdf;

                             L += day_L * pr;
                          }
                          cum_e *= ee[r];           
                       }
                       if(L <= 0.0)
                       {
                          printf("community_log_likelihood (line 414):\n");
                          printf("Infection likelihood must be positive!  L=%f\n", L);
                          printf("Person id=%d  start_day=%d  day_ill=%d  inf1=%d  inf2=%d\n", 
                                  person->id, start_day, person->day_ill, inf1, inf2);
                          for(k=0; k<n_b_mode; k++)  printf("b[%d]=%e  ", k, b[k]);   
                          for(k=0; k<n_p_mode; k++)  printf("p[%d]=%e  ", k, p[k]);
                          printf("\n");     
                          for(t=inf1;t<=inf2;t++)
                          {
                             r = t - start_day;
                             printf("ee[]=%e 1-ee[]=%e\n", ee[r], 1.0-ee[r]);
                          } 
                          error = 1;
                          exit(0);
                       }
                       log_L = log(L);
                       r = inf1 - start_day;
                       if(r > 0)  log_L += cum_log_ee[r-1];
                    }
                 }
                 /******************************************************************************************
                  calculate likelihood for an escaped subject, but need to consider right censoring. 
                  ******************************************************************************************/
                 else  /* if never been infected */
                 {
                    if(cfg_pars.adjust_for_right_censoring == 0)
                    {   
                       r = stop_day - start_day;
                       log_L = cum_log_ee[r];
                    }
                    else
                    {
                       inf1 = stop_day - cfg_pars.max_incubation + 1;
                       inf2 = stop_day - cfg_pars.min_incubation;
                       if(inf2 >= start_day)
                       {
                          // it could happen that inf1 < start_day < inf2
                          inf1 = max(inf1, start_day);

                          if(inf2 >= inf1)
                          {
                             cum_e = 1.0;
                             //if(person->id == 1)  printf("inf1=%d  inf2=%d\n", inf1, inf2);
                             // potential infection during days from inf1 to inf2
                             for(t=inf1; t<=inf2; t++)
                             {
                                pr = sdf_incubation[inf2 - t];
                                r = t - start_day;
                                if(ee[r] < 0.0 || ee[r] > 1.0)
                                {
                                   printf("community_loglikelihood (line 464), escaped person id=%d: t=%d  ee[r]=%e\n", person->id, t, ee[r]);
                                   error = 1;
                                   exit(0);
                                }
                                if(ee[r] >= 0.0 && ee[r] < 1.0)
                                {
                                   p_inf = (fabs(log(ee[r]))<1e-6)? (-log(ee[r])) : (1.0 - ee[r]);
                                   day_L = cum_e * p_inf;
                                   L += day_L * pr;
                                }
                                cum_e *= ee[r];
                                //if(person->id == 1)  printf("t=%d  ee=%e  p_inf=%e  cum_e=%e  day_L=%e  pr=%e  L=%e\n", t, ee[r], p_inf, cum_e, day_L, pr, L);
                             }

                             // escape during days from inf1 to inf2 
                             day_L = cum_e;
                             L += day_L;
                             if(L <= 0.0)
                             {
                                printf("community_loglikelihood (line 483): escaped person id=%d, L=%e is not positive\n", person->id, L);
                                error = 1;
                                exit(0);
                             }

                             log_L = log(L);
                          }   
                          r = inf1 - start_day;
                          if(r > 0) log_L += cum_log_ee[r-1];
                       }
                    }
                 } 
                 my_log_L += log_L;
              }
           } //if( person->pre_immune == 0)
           //printf("person %d: my_log_L=%e\n", person->id, my_log_L);
           log_L_all += person->weight * my_log_L;
           ptr_integer = ptr_integer->next;
        } //while(ptr_integer != NULL)
        ptr_class = ptr_class->next;
     } //while(ptr_class != NULL)
  }

end:
  if(error == 1)  log_L_all = -1e200;
  free(par); 

  return(log_L_all);
}

/* this program outputs log-likelihood and first derivatives for a given community
and given parameters. It is mainly used for checking MC error*/
int community_score(COMMUNITY *community, double *par_effective, double *log_likelihood, MATRIX *first)
{
  int i, j, h, k, l, m, n;
  int representative, converge, error=0;
  int start_day, stop_day, day_last_obs, min_size;
  int inf1, inf2;
  int r, rr, t, id, index, skip, verify;
  int found, b_mode, p_mode, q_mode, u_mode;
  int n_par, n_b_mode, n_p_mode, n_q_mode, n_u_mode;
  int n_c2p_covariate, n_p2p_covariate, n_imm_covariate, n_pat_covariate;
  int n_covariate, n_time_ind_covariate, n_time_dep_covariate;
  int n_sus_p2p_covariate, n_inf_p2p_covariate, n_int_p2p_covariate; 
  int n_par_equiclass;
  int stop, positive, n_iter;
  int shift, shift_row, shift_col;

  double factor, factor1, factor2, factor3;
  double e, f, ff, s, p_inf, Q, U;
  double covariate_effect, cum_e, day_L, L, LQ, log_L, my_log_L, my_L, log_L_all;
  
  double *par, sdf, pdf, pr;
  double logit_f;
  double temp, asym_effect;

  RISK *ptr_risk;
  PEOPLE *person, *member;
  RISK_CLASS *ptr_class;
  INTEGER_CHAIN *ptr_integer;
  

  /***********************************************************************************************************
   * n_b_mode: number of types/modes of community-to-person contact. 
   * n_p_mode: number of types/modes of person-to-person contact.
   * n_time_ind_covariate: number of time-independent covariates
   * n_time_dep_covariate: number of time-dependent covariates
   * n_covariate: total number of covariates
   * n_c2p_covariate: number of covariates of the susceptible person that can modify community-to-person risk
   * n_sus_p2p_covariate: number of covariates of the susceptible person that can modify person-to-person risk
   * n_inf_p2p_covariate: number of covariates of the infective person that can modify person-to-person risk
   * n_inf_p2p_covariate: number of interactions between covariates of the susceptible person
   *                      and the infective person that can modify person-to-person risk
   * n_p2p_covariate: total number of covariates that can modify person-to-person risk, which is the sum
   *                  of above three.
   * The user need to supply  n_time_ind_covariate, n_time_dep_covariate, n_c2p_covariate,
   * n_sus_p2p_covariate, n_inf_p2p_covariate and n_int_p2p_covariate.
   ***********************************************************************************************************/
  n_b_mode = cfg_pars.n_b_mode;
  n_p_mode = cfg_pars.n_p_mode;
  n_u_mode = cfg_pars.n_u_mode;
  n_q_mode = cfg_pars.n_q_mode;
  n_time_ind_covariate = cfg_pars.n_time_ind_covariate;
  n_time_dep_covariate = cfg_pars.n_time_dep_covariate;
  n_covariate = cfg_pars.n_covariate;
  n_c2p_covariate = cfg_pars.n_c2p_covariate;
  n_sus_p2p_covariate = cfg_pars.n_sus_p2p_covariate;
  n_inf_p2p_covariate = cfg_pars.n_inf_p2p_covariate;
  n_int_p2p_covariate = cfg_pars.n_int_p2p_covariate;
  n_p2p_covariate = cfg_pars.n_p2p_covariate;
  n_pat_covariate = cfg_pars.n_pat_covariate;
  n_imm_covariate = cfg_pars.n_imm_covariate;
  n_par = cfg_pars.n_par;
  n_par_equiclass = cfg_pars.n_par_equiclass;
  asym_effect = cfg_pars.asym_effect_est;
  

  make_1d_array_double(&par, n_par, 0.0);


  /************************************************************************************************************
   * f is person-to-person escape probability. e is the daily escape probability for a susceptible.
   * log_f_lb is the 1st derivative of log(f) with respect to lb.
   * log_f_lb_lb is the second derivative of log(f) with respect to lb.
   * Other derivative terms are similarly defined.
   *
   * e is the daily escape probability for a susceptible from all contacts.
   * cum_e and cum_log_e_* are cumulatives of e and log_e_* over a short period, 
   * from day_ill-max_latent to t, where day_ill-max_latent<=t<=day_ill-min_latent.
   * 
   * ee[t] is the escape probability for day t. 
   * cum_log_ee[t] = log(ee[1] * ee[2] * ... * ee[t])
   * 
   * temp_* variables are just auxiliary variables to help get cum_log_ee_* variables.
   * 
   * day_L is the likelihood for each infected person for the period 
   * from day_ill-max_latent+1 to t, where day_ill-max_latent<=t<=day_ill-min_latent.
   * log_day_L_* are derivatives of log(day_L). 
   *
   * L is the overall likelihood, and L_* and log_L_* are derivatives of L and log(L) respectively.
   * *********************************************************************************************************/

  //make sure input parameter values are not out of scope
  for(i=0; i<n_par_equiclass; i++)
  {
     if(par_effective[i] < cfg_pars.lower_search_bound[i])
     {
        par_effective[i] = cfg_pars.lower_search_bound[i];
     }
     if(par_effective[i] > cfg_pars.upper_search_bound[i])
     {
        par_effective[i] = cfg_pars.upper_search_bound[i];
     }
     m = cfg_pars.par_equiclass[i].member[0] - 1;
     if(m < n_b_mode + n_p_mode)
     {
        if(par_effective[i] < logit(close_to_0))  par_effective[i] = logit(close_to_0);
        if(par_effective[i] > logit(close_to_05))  par_effective[i] = logit(close_to_05);
     }
     else if(m < n_b_mode + n_p_mode + n_u_mode + n_q_mode)
     {
        if(par_effective[i] < logit(close_to_0))  par_effective[i] = logit(close_to_0);
        if(par_effective[i] > logit(close_to_1))  par_effective[i] = logit(close_to_1);
     }
  }

  //printf("Inside community_score()\n");
  //show_par_effective(par_effective, NULL);

  for(i=0; i<n_par_equiclass; i++)
  {
     for(j=0; j<cfg_pars.par_equiclass[i].size; j++)
     {
        m = cfg_pars.par_equiclass[i].member[j] - 1;
        par[m] = par_effective[i];
     }
  }

  //if fixed parameters exist, need to set fixed parameters to prespecified values
  if(cfg_pars.n_par_fixed > 0)
  {
     for(i=0; i<cfg_pars.n_par_fixed; i++)
     {
        j = cfg_pars.par_fixed_id[i] - 1;
        par[j] =  cfg_pars.par_fixed_value[i];
     }
  }
  
  
  for(i=0; i<n_b_mode; i++)  {lb[i] = par[i]; b[i] = inv_logit(lb[i]);}  
  for(i=0; i<n_p_mode; i++)  {lp[i] = par[n_b_mode + i]; p[i] = inv_logit(lp[i]);}  
  for(i=0; i<n_u_mode; i++)  {lu[i] = par[n_b_mode + n_p_mode + i]; u[i] = inv_logit(lu[i]);}  
  for(i=0; i<n_q_mode; i++)  {lq[i] = par[n_b_mode + n_p_mode + n_u_mode + i]; q[i] = inv_logit(lq[i]);}  
  for(i=0; i<n_c2p_covariate; i++)  coeff_c2p[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode + i];  
  for(i=0; i<n_p2p_covariate; i++)  coeff_p2p[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode 
                                                 + n_c2p_covariate + i]; 
  for(i=0; i<n_pat_covariate; i++)  coeff_pat[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode 
                                                 + n_c2p_covariate + n_p2p_covariate + i]; 
  for(i=0; i<n_imm_covariate; i++)  coeff_imm[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode 
                                                 + n_c2p_covariate + n_p2p_covariate + n_pat_covariate + i]; 
  
  log_L_all = 0.0;
  for(j=0; j<n_b_mode; j++) score_lb[j] = 0.0;
  for(j=0; j<n_p_mode; j++) score_lp[j] = 0.0;
  for(j=0; j<n_u_mode; j++) score_lu[j] = 0.0;
  for(j=0; j<n_q_mode; j++) score_lq[j] = 0.0;
  for(k=0; k<n_c2p_covariate; k++) score_c2p[k] = 0.0;
  for(k=0; k<n_p2p_covariate; k++) score_p2p[k] = 0.0;
  for(k=0; k<n_pat_covariate; k++) score_pat[k] = 0.0;
  for(k=0; k<n_imm_covariate; k++) score_imm[k] = 0.0;

  if(cfg_pars.adjust_for_left_truncation == 1)  min_size = community->size - community->size_idx;
  else  min_size = community->size;
  if(min_size > 0 && community->ignore == 0)
  {
     // If case-ascertained design, determine the earliest and latest index illness onset of index cases
     // which is used to adjust for left truncation (selection bias).
     // The following code is disabled because it has been done in set_index
     /*if(cfg_pars.adjust_for_left_truncation == 1)
     {
        i = community->idx[0];
        community->earliest_idx_day_ill = community->latest_idx_day_ill = people[i].day_ill;
        for(j=1; j<community->size_idx; j++)
        {
           i = community->idx[j];
           if(people[i].day_ill < community->earliest_idx_day_ill)
              community->earliest_idx_day_ill = people[i].day_ill;
           if(people[i].day_ill > community->latest_idx_day_ill)
              community->latest_idx_day_ill = people[i].day_ill;
        }
     }*/

     // index cases are now included in the risk classes.
     ptr_class = community->risk_class;
     while(ptr_class != NULL)
     {

        if(cfg_pars.common_contact_history_within_community == 1)
        {
           temp = 0.0;
           for(j=0; j<n_b_mode; j++) temp_lb[j] = 0.0;
           for(j=0; j<n_p_mode; j++) temp_lp[j] = 0.0;
           for(k=0; k<n_c2p_covariate; k++) temp_c2p[k] = 0.0;
           for(k=0; k<n_p2p_covariate; k++) temp_p2p[k] = 0.0;

           if(cfg_pars.adjust_for_left_truncation == 1 && community->earliest_idx_day_ill != MISSING)
              start_day = max(community->day_epi_start, community->earliest_idx_day_ill - cfg_pars.max_incubation + 1);
           else  start_day = community->day_epi_start;
           stop_day = community->day_epi_stop;
           //printf("community=%d  day_epi_start=%d earliest_idx_day_ill=%d \n", community->id, community->day_epi_start, community->earliest_idx_day_ill);
           //printf("id=%d start_day=%d  stop_day=%d\n", person->id, start_day, stop_day);
           for(t=start_day; t<=stop_day; t++)
           {
              r = t - start_day;
              rr = t - community->day_epi_start; //the time reference for risk history is community->day_epi_start, not start_day.
              e = 1.0;
              for(j=0; j<n_b_mode; j++) log_e_lb[j] = 0.0;
              for(j=0; j<n_p_mode; j++) log_e_lp[j] = 0.0;
              for(k=0; k<n_c2p_covariate; k++) log_e_c2p[k] = 0.0;
              for(k=0; k<n_p2p_covariate; k++) log_e_p2p[k] = 0.0;

              //calculate log_f_* variables related to the probability of escaping risk from common source.
              //update log_e_* variables accordingly
              ptr_risk = ptr_class->risk_history[rr].c2p_risk;
              while(ptr_risk != NULL)
              {
                 covariate_effect = 0.0;
                 for(k=0; k<n_c2p_covariate; k++)
                    covariate_effect += ptr_risk->covariate[k] * coeff_c2p[k];
                 covariate_effect += ptr_risk->offset;
                 
                 b_mode = ptr_risk->contact_mode;
                 logit_f = lb[b_mode] + covariate_effect;
                 
                 ff = inv_logit(logit_f);
                 f = 1.0 - ff;
                 log_f_lb[b_mode] = - ff;
                 
                 for(k=0; k<n_c2p_covariate; k++)
                    log_f_c2p[k] = -ptr_risk->covariate[k] * ff;
                 
                 e *= ipow(f, ptr_risk->size);
                 log_e_lb[b_mode] += ptr_risk->size * log_f_lb[b_mode];
                 for(k=0; k<n_c2p_covariate; k++)
                    log_e_c2p[k] += ptr_risk->size * log_f_c2p[k];
                 ptr_risk = ptr_risk->next;
              }

              //calculate log_f_* variables related to the probability of escaping risk from infective people.
              //update log_e_* variables accordingly
              //The difference in infectiousness level between symptomatic and asymptomatic cases
              //is adjusted by one additional covariate that denote the symptom status of the infective person.
              ptr_risk = ptr_class->risk_history[rr].p2p_risk;
              while(ptr_risk != NULL)
              {
                 covariate_effect = 0.0;
                 for(k=0; k<n_p2p_covariate; k++)
                    covariate_effect += ptr_risk->covariate[k] * coeff_p2p[k];
                 covariate_effect += ptr_risk->offset;
                 
                 p_mode = ptr_risk->contact_mode;
                 logit_f = lp[p_mode] + covariate_effect;
                 s = ptr_risk->infective_prob;
                 if( s > 0)
                 {
                    ff = inv_logit(logit_f) * s * bipow(asym_effect, 1 - ptr_risk->symptom);
                    f = 1 - ff;
                    factor1 = ff / s;
                    log_f_lp[p_mode] = - ff * (1 - factor1) / f;
                    
                    for(k=0; k<n_p2p_covariate; k++)\
                       log_f_p2p[k] = log_f_lp[p_mode] * ptr_risk->covariate[k];
                    
                    if(f <= 0.0 || f > 1.0)
                    {
                       printf("community_score (line 776): f=%e\n", f);
                       error = 1;
                       goto end;
                    }
                    e *= ipow(f, ptr_risk->size);
                    log_e_lp[p_mode] += ptr_risk->size * log_f_lp[p_mode];
                    for(k=0; k<n_p2p_covariate; k++)
                       log_e_p2p[k] += ptr_risk->size * log_f_p2p[k];
                 }
                 ptr_risk = ptr_risk->next;
              }

              //update log_ee_* variables, which are related the cumulative escape probabilities
              ee[r] = e;
              for(j=0; j<n_b_mode; j++) log_ee_lb[r][j] = log_e_lb[j];
              for(j=0; j<n_p_mode; j++) log_ee_lp[r][j] = log_e_lp[j];
              for(k=0; k<n_c2p_covariate; k++) log_ee_c2p[r][k] = log_e_c2p[k];
              for(k=0; k<n_p2p_covariate; k++) log_ee_p2p[r][k] = log_e_p2p[k];

              //printf("t=%d  log_e_lb=%e\n", t, log_e_lb[0]);
              temp += log(ee[r]);
              cum_log_ee[r] = temp;

              for(j=0; j<n_b_mode; j++)
              {
                 temp_lb[j] += log_ee_lb[r][j];
                 cum_log_ee_lb[r][j] = temp_lb[j];
              }

              for(j=0; j<n_p_mode; j++)
              {
                 temp_lp[j] += log_ee_lp[r][j];
                 cum_log_ee_lp[r][j] = temp_lp[j];
              }

              for(k=0; k<n_c2p_covariate; k++)
              {
                 temp_c2p[k] += log_ee_c2p[r][k];
                 cum_log_ee_c2p[r][k] = temp_c2p[k];
              }   

              for(k=0; k<n_p2p_covariate; k++)
              {
                 temp_p2p[k] += log_ee_p2p[r][k];
                 cum_log_ee_p2p[r][k] = temp_p2p[k];
              }   
           } /* end of day t */
        }// if(cfg_pars.common_contact_history_within_community == 1)

        ptr_integer = ptr_class->member;
        while(ptr_integer != NULL)
        {
           person = people + ptr_integer->id;
           if(cfg_pars.adjust_for_left_truncation == 1 &&  person->idx == 1)  
           {
              ptr_integer = ptr_integer->next;
              continue;
           }   
           if(n_u_mode > 0)  u_mode = person->u_mode;
           if(n_q_mode > 0)  q_mode = person->q_mode;
           
           my_log_L = 0.0;
           for(j=0; j<n_b_mode; j++) my_log_L_lb[j] = 0;
           for(j=0; j<n_p_mode; j++) my_log_L_lp[j] = 0;
           for(k=0; k<n_c2p_covariate; k++) my_log_L_c2p[k] = 0;
           for(k=0; k<n_p2p_covariate; k++) my_log_L_p2p[k] = 0;

           if(n_u_mode > 0)
           {
              my_log_L_lu = 0;
              
              for(k=0; k<n_pat_covariate; k++)
                 my_log_L_pat[k] = 0;
           }

           if(n_q_mode > 0)
           {
              my_log_L_lq = 0;
              
              for(k=0; k<n_imm_covariate; k++)
                 my_log_L_imm[k] = 0;
           }

           if(n_q_mode > 0)
           {
              covariate_effect = 0.0;
              for(k=0; k<n_imm_covariate; k++)
                 covariate_effect += person->imm_covariate[k] * coeff_imm[k];
              Q = inv_logit(lq[q_mode] + covariate_effect);
              Q_lq = Q * (1 - Q);
              for(k=0; k<n_imm_covariate; k++)
                 Q_imm[k] = person->imm_covariate[k] * Q * (1 - Q);

              if( person->pre_immune == 1)
              {
                 // consider prior immunity levels
                    factor1 = 1 / Q;
                    factor2 = -factor1 * factor1;
                    my_log_L += log(Q);
              }
              else
              {
                 factor1 = -1 / (1 - Q);
                 factor2 = factor1 / (1 - Q);
                 my_log_L += log(1 - Q);
              }
              
              my_log_L_lq += factor1 * Q_lq;

              for(k=0; k<n_imm_covariate; k++)
                 my_log_L_imm[k] += factor1 * Q_imm[k];
           }


           // if people in the comunity do not share the same exposure/risk history,
           // calculate individual-level exposure/risk history.
           if(cfg_pars.common_contact_history_within_community == 0 && person->pre_immune == 0)
           {
              temp = 0.0;
              for(j=0; j<n_b_mode; j++) temp_lb[j] = 0.0;
              for(j=0; j<n_p_mode; j++) temp_lp[j] = 0.0;
              for(k=0; k<n_c2p_covariate; k++) temp_c2p[k] = 0.0;
              for(k=0; k<n_p2p_covariate; k++) temp_p2p[k] = 0.0;

              if(cfg_pars.adjust_for_left_truncation == 1 && community->earliest_idx_day_ill != MISSING)
                 start_day = max(community->day_epi_start, community->earliest_idx_day_ill - cfg_pars.max_incubation + 1);
              else  start_day = community->day_epi_start;
              day_last_obs = min(community->day_epi_stop, person->day_exit);
              day_last_obs = min(community->day_last_followup, day_last_obs);
              stop_day = (person->infection == 1)? person->day_infection_upper:day_last_obs;
              if(stop_day < start_day)
              {
                 printf("community_score (line 908): community %d, person %d, stop day %d < start day %d\n", community->id, person->id, stop_day, start_day);
                 exit(0);
              }
              for(t=start_day; t<=stop_day; t++)
              {
                 r = t - start_day;
                 rr = t - community->day_epi_start; //the time reference for risk history is community->day_epi_start, not start_day.
                 e = 1.0;
                 for(j=0; j<n_b_mode; j++) log_e_lb[j] = 0.0;
                 for(j=0; j<n_p_mode; j++) log_e_lp[j] = 0.0;
                 for(k=0; k<n_c2p_covariate; k++) log_e_c2p[k] = 0.0;
                 for(k=0; k<n_p2p_covariate; k++) log_e_p2p[k] = 0.0;

                 //calculate log_f_* variables related to the probability of escaping risk from common source.
                 //update log_e_* variables accordingly
                 ptr_risk = person->risk_history[rr].c2p_risk;
                 while(ptr_risk != NULL)
                 {
                    covariate_effect = 0.0;
                    for(k=0; k<n_c2p_covariate; k++)
                       covariate_effect += ptr_risk->covariate[k] * coeff_c2p[k];
                    covariate_effect += ptr_risk->offset;
                    
                    b_mode = ptr_risk->contact_mode;
                    logit_f = lb[b_mode] + covariate_effect;
                    
                    ff = inv_logit(logit_f);
                    f = 1.0 - ff;
                    log_f_lb[b_mode] = - ff;
                    
                    for(k=0; k<n_c2p_covariate; k++)
                       log_f_c2p[k] = -ptr_risk->covariate[k] * ff;
                    
                    e *= ipow(f, ptr_risk->size);
                    log_e_lb[b_mode] += ptr_risk->size * log_f_lb[b_mode];
                    for(k=0; k<n_c2p_covariate; k++)
                       log_e_c2p[k] += ptr_risk->size * log_f_c2p[k];
                    ptr_risk = ptr_risk->next;
                 }

                 //calculate log_f_* variables related to the probability of escaping risk from infective people.
                 //update log_e_* variables accordingly
                 //The difference in infectiousness level between symptomatic and asymptomatic cases
                 //is adjusted by one additional covariate that denote the symptom status of the infective person.
                 ptr_risk = person->risk_history[rr].p2p_risk;
                 while(ptr_risk != NULL)
                 {
                    covariate_effect = 0.0;
                    for(k=0; k<n_p2p_covariate; k++)
                       covariate_effect += ptr_risk->covariate[k] * coeff_p2p[k];
                    covariate_effect += ptr_risk->offset;
                    
                    p_mode = ptr_risk->contact_mode;
                    logit_f = lp[p_mode] + covariate_effect;
                    s = ptr_risk->infective_prob;
                    if( s > 0)
                    {
                       ff = inv_logit(logit_f) * s * bipow(asym_effect, 1 - ptr_risk->symptom);
                       f = 1 - ff;
                       factor1 = ff / s;
                       log_f_lp[p_mode] = - ff * (1 - factor1) / f;
                       
                       for(k=0; k<n_p2p_covariate; k++)
                          log_f_p2p[k] = log_f_lp[p_mode] * ptr_risk->covariate[k];
                       
                       if(f <= 0.0 || f > 1.0)
                       {
                          printf("Community_score (line 975): f=%e\n", f);
                          error = 1;
                          goto end;
                       }
                       e *= ipow(f, ptr_risk->size);
                       log_e_lp[p_mode] += ptr_risk->size * log_f_lp[p_mode];
                       for(k=0; k<n_p2p_covariate; k++)
                          log_e_p2p[k] += ptr_risk->size * log_f_p2p[k];
                    }
                    ptr_risk = ptr_risk->next;
                 }

                 //update log_ee_* variables, which are related the cumulative escape probabilities
                 ee[r] = e;
                 for(j=0; j<n_b_mode; j++) log_ee_lb[r][j] = log_e_lb[j];
                 for(j=0; j<n_p_mode; j++) log_ee_lp[r][j] = log_e_lp[j];
                 for(k=0; k<n_c2p_covariate; k++) log_ee_c2p[r][k] = log_e_c2p[k];
                 for(k=0; k<n_p2p_covariate; k++) log_ee_p2p[r][k] = log_e_p2p[k];


                 temp += log(ee[r]);
                 cum_log_ee[r] = temp;

                 for(j=0; j<n_b_mode; j++)
                 {
                    temp_lb[j] += log_ee_lb[r][j];
                    cum_log_ee_lb[r][j] = temp_lb[j];
                 }

                 for(j=0; j<n_p_mode; j++)
                 {
                    temp_lp[j] += log_ee_lp[r][j];
                    cum_log_ee_lp[r][j] = temp_lp[j];
                 }

                 for(k=0; k<n_c2p_covariate; k++)
                 {
                    temp_c2p[k] += log_ee_c2p[r][k];
                    cum_log_ee_c2p[r][k] = temp_c2p[k];
                 }   

                 for(k=0; k<n_p2p_covariate; k++)
                 {
                    temp_p2p[k] += log_ee_p2p[r][k];
                    cum_log_ee_p2p[r][k] = temp_p2p[k];
                 }   
              } /* end of day t */
           }

           if( person->pre_immune == 0)
           {
              if(cfg_pars.adjust_for_left_truncation == 1 && community->earliest_idx_day_ill != MISSING)
                 start_day = max(community->day_epi_start, community->earliest_idx_day_ill - cfg_pars.max_incubation + 1);
              else  start_day = community->day_epi_start;

              day_last_obs = min(community->day_epi_stop, person->day_exit);
              day_last_obs = min(community->day_last_followup, day_last_obs);
              stop_day = (person->infection == 1)? person->day_infection_upper:day_last_obs;

              if(stop_day < start_day)
              {
                 printf("community_score (line 1036): community %d, person %d, stop day %d < start day %d\n", community->id, person->id, stop_day, start_day);
                 exit(0);
              }
              else
              {
                 L = 0.0;
                 for(j=0; j<n_b_mode; j++) L_lb[j] = 0.0;
                 for(j=0; j<n_p_mode; j++) L_lp[j] = 0.0;
                 for(k=0; k<n_c2p_covariate; k++) L_c2p[k] = 0.0;
                 for(k=0; k<n_p2p_covariate; k++) L_p2p[k] = 0.0;
                 
                 if(n_u_mode > 0)
                 {
                    L_lu = 0.0;
                    for(k=0; k<n_pat_covariate; k++)
                       L_pat[k] = 0.0;
                 }

                 log_L = 0.0;
                 for(j=0; j<n_b_mode; j++) log_L_lb[j] = 0.0;
                 for(j=0; j<n_p_mode; j++) log_L_lp[j] = 0.0;
                 for(k=0; k<n_c2p_covariate; k++) log_L_c2p[k] = 0.0;
                 for(k=0; k<n_p2p_covariate; k++) log_L_p2p[k] = 0.0;
                 
                 if(n_u_mode > 0)
                 {
                    log_L_lu = 0.0;
                    for(k=0; k<n_pat_covariate; k++)
                       log_L_pat[k] = 0.0;
                 }


                 /******************************************************************************************
                  calculate likelihood for a infected case with symptoms 
                  ******************************************************************************************/
                    
                 if(person->infection == 1) /* if he was infected */
                 {
                    inf1 = person->day_infection_lower;
                    inf2 = person->day_infection_upper;
                    if(inf2 >= start_day)
                    {
                       // it could happen that inf1 < start_day < inf2
                       inf1 = max(inf1, start_day);

                       cum_e = 1.0; 
                       for(j=0; j<n_b_mode; j++) cum_log_e_lb[j] = 0.0; 
                       for(j=0; j<n_p_mode; j++) cum_log_e_lp[j] = 0.0; 
                       for(k=0; k<n_c2p_covariate; k++) cum_log_e_c2p[k] = 0.0; 
                       for(k=0; k<n_p2p_covariate; k++) cum_log_e_p2p[k] = 0.0; 

                       for(t=inf1; t<=inf2; t++)
                       {
                          pdf = pdf_incubation[person->day_ill - cfg_pars.min_incubation - t];
                          r = t - start_day;
                          if(ee[r] < 0.0 || ee[r] > 1.0)
                          {
                             printf("community_score (line 1093): infected person id=%d: t=%d ee[r]=%e\n", person->id, t, ee[r]);
                             error = 1;
                             goto end;
                          }
                             
                          if(ee[r] >= 0.0 && ee[r] < 1.0)
                          {
                             p_inf = (fabs(log(ee[r]))<1e-6)? (-log(ee[r])) : (1.0 - ee[r]);
                             factor1 = -ee[r] / p_inf; 
                             day_L = cum_e * p_inf; 
                             for(j=0; j<n_b_mode; j++) 
                                log_day_L_lb[j] = cum_log_e_lb[j] + factor1 * log_ee_lb[r][j]; 
                             for(j=0; j<n_p_mode; j++) 
                                log_day_L_lp[j] = cum_log_e_lp[j] + factor1 * log_ee_lp[r][j]; 
                             for(k=0; k<n_c2p_covariate; k++) 
                                log_day_L_c2p[k] = cum_log_e_c2p[k] + factor1 * log_ee_c2p[r][k]; 
                             for(k=0; k<n_p2p_covariate; k++) 
                                log_day_L_p2p[k] = cum_log_e_p2p[k] + factor1 * log_ee_p2p[r][k]; 

                             // now consider pathogenicity
                             if(n_u_mode > 0)
                             {
                                rr = t - community->day_epi_start; // track covariate for pathogenicity, for which timeline starts from community->day_epi_start.
                                                                   // Covariates affecting pathogenicity are from the day of infection 
                                covariate_effect = 0.0;
                                for(k=0; k<n_pat_covariate; k++)
                                   covariate_effect += person->pat_covariate[k][rr] * coeff_pat[k];
                                U = inv_logit(lu[u_mode] + covariate_effect);
                                U_lu = U * (1 - U);

                                for(k=0; k<n_pat_covariate; k++)
                                   U_pat[k] = person->pat_covariate[k][rr] * U * (1 - U);
                                
                                pr = pdf * (person->symptom * U + (1 - person->symptom) * (1 - U));
                                pr_lu = bipow(-1, (1-person->symptom)) * pdf * U_lu;
                                for(k=0; k<n_pat_covariate; k++)  
                                   pr_pat[k] = bipow(-1, (1-person->symptom)) * pdf * U_pat[k];
                             }
                             else  pr = pdf;
         
                             temp = day_L * pr; 
                             L += temp; 
                             for(j=0; j<n_b_mode; j++) L_lb[j] += log_day_L_lb[j] * temp; 
                             for(j=0; j<n_p_mode; j++) L_lp[j] += log_day_L_lp[j] * temp; 
                             for(k=0; k<n_c2p_covariate; k++) L_c2p[k] += log_day_L_c2p[k] * temp; 
                             for(k=0; k<n_p2p_covariate; k++) L_p2p[k] += log_day_L_p2p[k] * temp; 
                             if(n_u_mode > 0)
                             {
                                L_lu += day_L * pr_lu;
                                for(k=0; k<n_pat_covariate; k++)
                                   L_pat[k] += day_L * pr_pat[k];
                             }
                          }

                          cum_e *= ee[r]; 
                          for(j=0; j<n_b_mode; j++) cum_log_e_lb[j] += log_ee_lb[r][j]; 
                          for(j=0; j<n_p_mode; j++) cum_log_e_lp[j] += log_ee_lp[r][j]; 
                          for(k=0; k<n_c2p_covariate; k++) cum_log_e_c2p[k] += log_ee_c2p[r][k]; 
                          for(k=0; k<n_p2p_covariate; k++) cum_log_e_p2p[k] += log_ee_p2p[r][k]; 
                       }
                       if(L <= 0.0)
                       {
                          printf("community_score (line 1155):\n");
                          printf("Infection likelihood must be positive!  L=%f\n", L);
                          printf("Person id=%d  start_day=%d  day_ill=%d  inf1=%d  inf2=%d\n", 
                                  person->id, start_day, person->day_ill, inf1, inf2);
                          for(k=0; k<n_b_mode; k++)  printf("b[%d]=%e  ", k, b[k]);   
                          for(k=0; k<n_p_mode; k++)  printf("p[%d]=%e  ", k, p[k]);
                          printf("\n");     
                          for(t=inf1;t<=inf2;t++)
                          {
                             r = t - start_day;
                             printf("ee[]=%e 1-ee[]=%e\n", ee[r], 1.0-ee[r]);
                          } 
                            
                          error = 1;
                          goto end;
                       }

                       factor1 = 1.0 / L; 
                       log_L = log(L); 
                       for(j=0; j<n_b_mode; j++) log_L_lb[j] = L_lb[j] * factor1; 
                       for(j=0; j<n_p_mode; j++) log_L_lp[j] = L_lp[j] * factor1; 
                       for(k=0; k<n_c2p_covariate; k++) log_L_c2p[k] = L_c2p[k] * factor1; 
                       for(k=0; k<n_p2p_covariate; k++) log_L_p2p[k] = L_p2p[k] * factor1; 

                       if(n_u_mode > 0)
                       {
                          log_L_lu = L_lu * factor1;
                          for(k=0; k<n_pat_covariate; k++)
                             log_L_pat[k] = L_pat[k] * factor1;
                       }

                       r = inf1 - start_day; 
                       if(r > 0) 
                       { 
                          log_L += cum_log_ee[r-1]; 
                          for(j=0; j<n_b_mode; j++) 
                             log_L_lb[j] += cum_log_ee_lb[r-1][j]; 
                          for(j=0; j<n_p_mode; j++) 
                             log_L_lp[j] += cum_log_ee_lp[r-1][j]; 
                          for(k=0; k<n_c2p_covariate; k++) 
                             log_L_c2p[k] += cum_log_ee_c2p[r-1][k]; 
                          for(k=0; k<n_p2p_covariate; k++) 
                             log_L_p2p[k] += cum_log_ee_p2p[r-1][k]; 
                       }
                    }
                 }
                 /******************************************************************************************
                  calculate likelihood for an escaped subject, but need to consider right censoring. 
                  ******************************************************************************************/
                 else  /* if never been infected */
                 {
                       
                    if(cfg_pars.adjust_for_right_censoring == 0)
                    {
                       r = stop_day - start_day;
                       log_L = cum_log_ee[r];
                       for(j=0; j<n_b_mode; j++)
                          log_L_lb[j] += cum_log_ee_lb[r][j];
                       for(j=0; j<n_p_mode; j++)
                          log_L_lp[j] += cum_log_ee_lp[r][j];
                       for(k=0; k<n_c2p_covariate; k++)
                          log_L_c2p[k] += cum_log_ee_c2p[r][k];
                       for(k=0; k<n_p2p_covariate; k++)
                          log_L_p2p[k] += cum_log_ee_p2p[r][k];
                    }
                    else
                    {
                       inf1 = stop_day - cfg_pars.max_incubation + 1;
                       inf2 = stop_day - cfg_pars.min_incubation;
                       if(inf2 >= start_day)
                       {
                          // it could happen that inf1 < start_day < inf2
                          inf1 = max(inf1, start_day);

                          if(inf2 >= inf1)
                          {
                             /* potential infection during days from inf1 to inf2*/
                             cum_e = 1.0; 
                             for(j=0; j<n_b_mode; j++) cum_log_e_lb[j] = 0.0; 
                             for(j=0; j<n_p_mode; j++) cum_log_e_lp[j] = 0.0; 
                             for(k=0; k<n_c2p_covariate; k++) cum_log_e_c2p[k] = 0.0; 
                             for(k=0; k<n_p2p_covariate; k++) cum_log_e_p2p[k] = 0.0; 

                             for(t=inf1; t<=inf2; t++)
                             {
                                pr = sdf_incubation[inf2 - t];
                                r = t - start_day;
                                if(ee[r] < 0.0 || ee[r] > 1.0)
                                {
                                   printf("community_score (line 1244): escapeed person id=%d: t=%d ee[r]=%e\n", person->id, t, ee[r]);
                                   error = 1;
                                   goto end;
                                }

                                if(ee[r] >= 0.0 && ee[r] < 1.0)
                                {
                                   p_inf = (fabs(log(ee[r]))<1e-6)? (-log(ee[r])) : (1.0 - ee[r]);
                                   factor1 = -ee[r] / p_inf; 
                                   day_L = cum_e * p_inf; 
                                   for(j=0; j<n_b_mode; j++) 
                                      log_day_L_lb[j] = cum_log_e_lb[j] + factor1 * log_ee_lb[r][j]; 
                                   for(j=0; j<n_p_mode; j++) 
                                      log_day_L_lp[j] = cum_log_e_lp[j] + factor1 * log_ee_lp[r][j]; 
                                   for(k=0; k<n_c2p_covariate; k++) 
                                      log_day_L_c2p[k] = cum_log_e_c2p[k] + factor1 * log_ee_c2p[r][k]; 
                                   for(k=0; k<n_p2p_covariate; k++) 
                                      log_day_L_p2p[k] = cum_log_e_p2p[k] + factor1 * log_ee_p2p[r][k]; 
                                   
                                   temp = day_L * pr; 
                                   L += temp; 
                                   for(j=0; j<n_b_mode; j++) 
                                      L_lb[j] += log_day_L_lb[j] * temp; 
                                   for(j=0; j<n_p_mode; j++) 
                                      L_lp[j] += log_day_L_lp[j] * temp; 
                                   for(k=0; k<n_c2p_covariate; k++) 
                                      L_c2p[k] += log_day_L_c2p[k] * temp; 
                                   for(k=0; k<n_p2p_covariate; k++) 
                                      L_p2p[k] += log_day_L_p2p[k] * temp; 
                                }
                                cum_e *= ee[r]; 
                                for(j=0; j<n_b_mode; j++) 
                                   cum_log_e_lb[j] += log_ee_lb[r][j]; 
                                for(j=0; j<n_p_mode; j++) 
                                   cum_log_e_lp[j] += log_ee_lp[r][j]; 
                                for(k=0; k<n_c2p_covariate; k++) 
                                   cum_log_e_c2p[k] += log_ee_c2p[r][k]; 
                                for(k=0; k<n_p2p_covariate; k++) 
                                   cum_log_e_p2p[k] += log_ee_p2p[r][k]; 
                             }

                             // escape during days from inf1 to inf2 
                             day_L = cum_e; 
                             for(j=0; j<n_b_mode; j++) 
                                log_day_L_lb[j] = cum_log_e_lb[j]; 
                             for(j=0; j<n_p_mode; j++) 
                                log_day_L_lp[j] = cum_log_e_lp[j]; 
                             for(k=0; k<n_c2p_covariate; k++) 
                                log_day_L_c2p[k] = cum_log_e_c2p[k]; 
                             for(k=0; k<n_p2p_covariate; k++) 
                                log_day_L_p2p[k] = cum_log_e_p2p[k]; 

                             L += day_L; 
                             for(j=0; j<n_b_mode; j++) 
                                L_lb[j] += log_day_L_lb[j] * day_L; 
                             for(j=0; j<n_p_mode; j++) 
                                L_lp[j] += log_day_L_lp[j] * day_L; 
                             for(k=0; k<n_c2p_covariate; k++) 
                                L_c2p[k] += log_day_L_c2p[k] * day_L; 
                             for(k=0; k<n_p2p_covariate; k++) 
                                L_p2p[k] += log_day_L_p2p[k] * day_L; 


                             if(L <= 0.0)
                             {
                                printf("community_score (line 1309): escaped person id=%d: L=%e is not positive\n", person->id, L);
                                error = 1;
                                goto end;
                             }
                             factor1 = 1.0 / L; 
                             log_L = log(L); 
                             for(j=0; j<n_b_mode; j++) 
                                log_L_lb[j] = L_lb[j] * factor1; 
                             for(j=0; j<n_p_mode; j++) 
                                log_L_lp[j] = L_lp[j] * factor1; 
                             for(k=0; k<n_c2p_covariate; k++) 
                                log_L_c2p[k] = L_c2p[k] * factor1; 
                             for(k=0; k<n_p2p_covariate; k++) 
                                log_L_p2p[k] = L_p2p[k] * factor1; 

                             if(n_u_mode > 0)
                             {
                                log_L_lu = L_lu * factor1;
                                for(k=0; k<n_pat_covariate; k++)
                                   log_L_pat[k] = L_pat[k] * factor1;
                             }
                          }
                          r = inf1 - start_day; 
                          if(r > 0) 
                          { 
                             log_L += cum_log_ee[r-1]; 
                             for(j=0; j<n_b_mode; j++) 
                                log_L_lb[j] += cum_log_ee_lb[r-1][j]; 
                             for(j=0; j<n_p_mode; j++) 
                                log_L_lp[j] += cum_log_ee_lp[r-1][j]; 
                             for(k=0; k<n_c2p_covariate; k++) 
                                log_L_c2p[k] += cum_log_ee_c2p[r-1][k]; 
                             for(k=0; k<n_p2p_covariate; k++) 
                                log_L_p2p[k] += cum_log_ee_p2p[r-1][k]; 
                          }
                       }
                    }                 
                 } 

                 my_log_L += log_L;
                 for(j=0; j<n_b_mode; j++)
                    my_log_L_lb[j] += log_L_lb[j];
                 for(j=0; j<n_p_mode; j++)
                    my_log_L_lp[j] += log_L_lp[j];
                 for(k=0; k<n_c2p_covariate; k++)
                    my_log_L_c2p[k] += log_L_c2p[k];
                 for(k=0; k<n_p2p_covariate; k++)
                    my_log_L_p2p[k] += log_L_p2p[k];

                 if(n_u_mode > 0)
                 {
                    my_log_L_lu += log_L_lu;
                    
                    for(k=0; k<n_pat_covariate; k++)
                       my_log_L_pat[k] += log_L_pat[k];
                 }
              }
           } //if(person->pre_immune == 0)

           for(j=0; j<n_b_mode; j++)
              score_lb[j] += person->weight * my_log_L_lb[j];
           for(j=0; j<n_p_mode; j++)
              score_lp[j] += person->weight * my_log_L_lp[j];
           for(j=0; j<n_c2p_covariate; j++)
              score_c2p[j] += person->weight * my_log_L_c2p[j];
           for(j=0; j<n_p2p_covariate; j++)
              score_p2p[j] += person->weight * my_log_L_p2p[j];

           if(n_u_mode > 0)
           {
              score_lu[u_mode] += person->weight * my_log_L_lu;
              for(k=0; k<n_pat_covariate; k++)
                 score_pat[k] += person->weight * my_log_L_pat[k];
           }      
              
           if(n_q_mode > 0)
           {
              score_lq[q_mode] += person->weight * my_log_L_lq;
              for(k=0; k<n_imm_covariate; k++)
                 score_imm[k] += person->weight * my_log_L_imm[k];
           }      

           log_L_all += person->weight * my_log_L;
           ptr_integer = ptr_integer->next;
        } //while(ptr_integer != NULL)
        ptr_class = ptr_class->next;
     } //while(ptr_class != NULL)
     //if(DEBUG == 1)  printf("%d: %e %e %e\n", person->id, score_lb[0], score_lp[0], score_lu[0]);
  } /*end if(community->size > 0)*/
  
  /* use logit transformation. newton-raphson algorithm is suitable for searching over the whole real line */
  shift = 0;
  for(i=0; i<n_b_mode; i++)  first->data[shift + i][0] = score_lb[i];
  shift = n_b_mode;
  for(i=0; i<n_p_mode; i++)  first->data[shift + i][0] = score_lp[i];
  shift = n_b_mode + n_p_mode;
  for(i=0; i<n_u_mode; i++)  first->data[shift + i][0] = score_lu[i];
  shift = n_b_mode + n_p_mode + n_u_mode;
  for(i=0; i<n_q_mode; i++)  first->data[shift + i][0] = score_lq[i];
  shift = n_b_mode + n_p_mode + n_u_mode + n_q_mode;
  for(i=0; i<n_c2p_covariate; i++)  first->data[shift + i][0] = score_c2p[i];
  shift = n_b_mode + n_p_mode + n_u_mode + n_q_mode + n_c2p_covariate;
  for(i=0; i<n_p2p_covariate; i++)  first->data[shift + i][0] = score_p2p[i];
  shift = n_b_mode + n_p_mode + n_u_mode + n_q_mode + n_c2p_covariate + n_p2p_covariate;
  for(i=0; i<n_pat_covariate; i++)  first->data[shift + i][0] = score_pat[i];
  shift = n_b_mode + n_p_mode + n_u_mode + n_q_mode + n_c2p_covariate + n_p2p_covariate + n_pat_covariate;
  for(i=0; i<n_imm_covariate; i++)  first->data[shift + i][0] = score_imm[i];

end:
  if(!(log_L_all <= 0))  log_L_all = -1e200;
  (* log_likelihood) = log_L_all;

  free(par); 
  return(error);
}


/* this program outputs log-likelihood, first and second derivatives for a given community
and given parameters. */


int community_derivatives(COMMUNITY *community, double *par_effective, double *log_likelihood, MATRIX *first, MATRIX *second)
{
  int i, j, h, k, l, m, n;
  int representative, converge, error=0;
  int start_day, stop_day, day_last_obs, min_size;
  int inf1, inf2;
  int r, rr, t, id, index, skip, verify;
  int found, b_mode, p_mode, q_mode, u_mode;
  int n_par, n_b_mode, n_p_mode, n_q_mode, n_u_mode;
  int n_c2p_covariate, n_p2p_covariate, n_imm_covariate, n_pat_covariate;
  int n_covariate, n_time_ind_covariate, n_time_dep_covariate;
  int n_sus_p2p_covariate, n_inf_p2p_covariate, n_int_p2p_covariate; 
  int n_par_equiclass;
  int stop, positive, n_iter;
  int shift, shift_row, shift_col;

  double factor, factor1, factor2, factor3;
  double e, f, ff, s, p_inf, Q, U;
  double covariate_effect, cum_e, day_L, L, LQ, log_L, my_log_L, my_L, log_L_all;
  
  double *par, sdf, pdf, pr;
  double logit_f;
  double temp, asym_effect;

  RISK *ptr_risk;
  PEOPLE *person, *member;
  RISK_CLASS *ptr_class;
  INTEGER_CHAIN *ptr_integer;
  

  /***********************************************************************************************************
   * n_b_mode: number of types/modes of community-to-person contact. 
   * n_p_mode: number of types/modes of person-to-person contact.
   * n_time_ind_covariate: number of time-independent covariates
   * n_time_dep_covariate: number of time-dependent covariates
   * n_covariate: total number of covariates
   * n_c2p_covariate: number of covariates of the susceptible person that can modify community-to-person risk
   * n_sus_p2p_covariate: number of covariates of the susceptible person that can modify person-to-person risk
   * n_inf_p2p_covariate: number of covariates of the infective person that can modify person-to-person risk
   * n_inf_p2p_covariate: number of interactions between covariates of the susceptible person
   *                      and the infective person that can modify person-to-person risk
   * n_p2p_covariate: total number of covariates that can modify person-to-person risk, which is the sum
   *                  of above three.
   * The user need to supply  n_time_ind_covariate, n_time_dep_covariate, n_c2p_covariate,
   * n_sus_p2p_covariate, n_inf_p2p_covariate and n_int_p2p_covariate.
   ***********************************************************************************************************/
  n_b_mode = cfg_pars.n_b_mode;
  n_p_mode = cfg_pars.n_p_mode;
  n_u_mode = cfg_pars.n_u_mode;
  n_q_mode = cfg_pars.n_q_mode;
  n_time_ind_covariate = cfg_pars.n_time_ind_covariate;
  n_time_dep_covariate = cfg_pars.n_time_dep_covariate;
  n_covariate = cfg_pars.n_covariate;
  n_c2p_covariate = cfg_pars.n_c2p_covariate;
  n_sus_p2p_covariate = cfg_pars.n_sus_p2p_covariate;
  n_inf_p2p_covariate = cfg_pars.n_inf_p2p_covariate;
  n_int_p2p_covariate = cfg_pars.n_int_p2p_covariate;
  n_p2p_covariate = cfg_pars.n_p2p_covariate;
  n_pat_covariate = cfg_pars.n_pat_covariate;
  n_imm_covariate = cfg_pars.n_imm_covariate;
  n_par = cfg_pars.n_par;
  n_par_equiclass = cfg_pars.n_par_equiclass;
  asym_effect = cfg_pars.asym_effect_est;
  

  make_1d_array_double(&par, n_par, 0.0);


  /************************************************************************************************************
   * f is person-to-person escape probability. e is the daily escape probability for a susceptible.
   * log_f_lb is the 1st derivative of log(f) with respect to lb.
   * log_f_lb_lb is the second derivative of log(f) with respect to lb.
   * Other derivative terms are similarly defined.
   *
   * e is the daily escape probability for a susceptible from all contacts.
   * cum_e and cum_log_e_* are cumulatives of e and log_e_* over a short period, 
   * from day_ill-max_latent to t, where day_ill-max_latent<=t<=day_ill-min_latent.
   * 
   * ee[t] is the escape probability for day t. 
   * cum_log_ee[t] = log(ee[1] * ee[2] * ... * ee[t])
   * 
   * temp_* variables are just auxiliary variables to help get cum_log_ee_* variables.
   * 
   * day_L is the likelihood for each infected person for the period 
   * from day_ill-max_latent+1 to t, where day_ill-max_latent<=t<=day_ill-min_latent.
   * log_day_L_* are derivatives of log(day_L). 
   *
   * L is the overall likelihood, and L_* and log_L_* are derivatives of L and log(L) respectively.
   * *********************************************************************************************************/

  //make sure input parameter values are not out of scope
  for(i=0; i<n_par_equiclass; i++)
  {
     if(par_effective[i] < cfg_pars.lower_search_bound[i])
     {
        par_effective[i] = cfg_pars.lower_search_bound[i];
     }
     if(par_effective[i] > cfg_pars.upper_search_bound[i])
     {
        par_effective[i] = cfg_pars.upper_search_bound[i];
     }
     m = cfg_pars.par_equiclass[i].member[0] - 1;
     if(m < n_b_mode + n_p_mode)
     {
        if(par_effective[i] < logit(close_to_0))  par_effective[i] = logit(close_to_0);
        if(par_effective[i] > logit(close_to_05))  par_effective[i] = logit(close_to_05);
     }
     else if(m < n_b_mode + n_p_mode + n_u_mode + n_q_mode)
     {
        if(par_effective[i] < logit(close_to_0))  par_effective[i] = logit(close_to_0);
        if(par_effective[i] > logit(close_to_1))  par_effective[i] = logit(close_to_1);
     }

  }

  for(i=0; i<n_par_equiclass; i++)
  {
     for(j=0; j<cfg_pars.par_equiclass[i].size; j++)
     {
        m = cfg_pars.par_equiclass[i].member[j] - 1;
        par[m] = par_effective[i];
     }
  }

  //if fixed parameters exist, need to set fixed parameters to prespecified values
  if(cfg_pars.n_par_fixed > 0)
  {
     for(i=0; i<cfg_pars.n_par_fixed; i++)
     {
        j = cfg_pars.par_fixed_id[i] - 1;
        par[j] =  cfg_pars.par_fixed_value[i];
     }
  }
  
  
  for(i=0; i<n_b_mode; i++)  {lb[i] = par[i]; b[i] = inv_logit(lb[i]);}  
  for(i=0; i<n_p_mode; i++)  {lp[i] = par[n_b_mode + i]; p[i] = inv_logit(lp[i]);}  
  for(i=0; i<n_u_mode; i++)  {lu[i] = par[n_b_mode + n_p_mode + i]; u[i] = inv_logit(lu[i]);}  
  for(i=0; i<n_q_mode; i++)  {lq[i] = par[n_b_mode + n_p_mode + n_u_mode + i]; q[i] = inv_logit(lq[i]);}  
  for(i=0; i<n_c2p_covariate; i++)  coeff_c2p[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode + i];  
  for(i=0; i<n_p2p_covariate; i++)  coeff_p2p[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode 
                                                 + n_c2p_covariate + i]; 
  for(i=0; i<n_pat_covariate; i++)  coeff_pat[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode 
                                                 + n_c2p_covariate + n_p2p_covariate + i]; 
  for(i=0; i<n_imm_covariate; i++)  coeff_imm[i] = par[n_b_mode + n_p_mode + n_u_mode + n_q_mode 
                                                 + n_c2p_covariate + n_p2p_covariate + n_pat_covariate + i]; 
  log_L_all = 0.0;
  INITIALIZE_SCORE_AND_INFO

  if(cfg_pars.adjust_for_left_truncation == 1)  min_size = community->size - community->size_idx;
  else  min_size = community->size;
  if(min_size > 0 && community->ignore == 0)
  {
     // If case-ascertained design, determine the earliest and latest index illness onset of index cases
     // which is used to adjust for left truncation (selection bias).
     // The following code is disabled because it has been done in set_index
     /*if(cfg_pars.adjust_for_left_truncation == 1)
     {
        i = community->idx[0];
        community->earliest_idx_day_ill = community->latest_idx_day_ill = people[i].day_ill;
        for(j=1; j<community->size_idx; j++)
        {
           i = community->idx[j];
           if(people[i].day_ill < community->earliest_idx_day_ill)
              community->earliest_idx_day_ill = people[i].day_ill;
           if(people[i].day_ill > community->latest_idx_day_ill)
              community->latest_idx_day_ill = people[i].day_ill;
        }
     }*/

     ptr_class = community->risk_class;
     while(ptr_class != NULL)
     {
        if(cfg_pars.common_contact_history_within_community == 1)
        {
           INITIALIZE_TEMP

           if(cfg_pars.adjust_for_left_truncation == 1 && community->earliest_idx_day_ill != MISSING)
              start_day = max(community->day_epi_start, community->earliest_idx_day_ill - cfg_pars.max_incubation + 1);
           else  start_day = community->day_epi_start;
           stop_day = community->day_epi_stop;
           if(stop_day < start_day)
           {
              printf("community_derivatives (line 1611): community %d, stop day %d < start day %d\n", community->id, stop_day, start_day);
              exit(0);
           }
           //printf("%d: start_day=%d  stop_day=%d\n", h, start_day, stop_day);
           for(t=start_day; t<=stop_day; t++)
           {
              r = t - start_day;
              rr = t - community->day_epi_start; //the time reference for risk history is community->day_epi_start, not start_day.
              INITIALIZE_LOG_E
              UPDATE_LOG_E_C2P
              //The difference in infectiousness level between symptomatic and asymptomatic cases
              //is adjusted by one additional covariate that denote the symptom status of the infective person.
              
              UPDATE_LOG_E_P2P
              
              UPDATE_LOG_EE
              
              UPDATE_CUM_LOG_EE

              //printf("t=%d  log_ee_lb=%e  log_ee_lb_lb=%e  log_ee_lp=%e  log_ee_lp_lp=%e\n", t, log_ee_lb[r][0], log_ee_lb_lb[r][0], log_ee_lp[r][0], log_ee_lp_lp[r][0]);
           } /* end of day t */
        }
        
        ptr_integer = ptr_class->member;
        while(ptr_integer != NULL)
        {
           person = people + ptr_integer->id;
           if(cfg_pars.adjust_for_left_truncation == 1 && person->idx == 1)  
           {
              ptr_integer = ptr_integer->next;
              continue;
           }   
     
           if(n_u_mode > 0)  u_mode = person->u_mode;
           if(n_q_mode > 0)  q_mode = person->q_mode;

           INITIALIZE_MY_LOG_L

           if(n_q_mode > 0)
           {
              covariate_effect = 0.0;
              for(k=0; k<n_imm_covariate; k++)
                 covariate_effect += person->imm_covariate[k] * coeff_imm[k];
              Q = inv_logit(lq[q_mode] + covariate_effect);
              Q_lq = Q * (1 - Q);
              Q_lq_lq = Q * (1 - Q) * (1 - 2 * Q);
              for(k=0; k<n_imm_covariate; k++)
                 Q_lq_imm[k] = person->imm_covariate[k] * Q * (1 - Q) * (1 - 2 * Q);
              for(k=0; k<n_imm_covariate; k++)
              {
                 Q_imm[k] = person->imm_covariate[k] * Q * (1 - Q);
                 for(l=k; l<n_imm_covariate; l++)
                    Q_imm_imm[k][l] = person->imm_covariate[k] * person->imm_covariate[l] * Q * (1 - Q) * (1 - 2 * Q);
              }

              if( person->pre_immune == 1)
              {
                 // consider prior immunity levels
                    factor1 = 1 / Q;
                    factor2 = -factor1 * factor1;
                    my_log_L += log(Q);
              }
              else
              {
                 factor1 = -1 / (1 - Q);
                 factor2 = factor1 / (1 - Q);
                 my_log_L += log(1 - Q);
              }
              
              my_log_L_lq += factor1 * Q_lq;
              my_log_L_lq_lq += factor1 * Q_lq_lq + factor2 * Q_lq * Q_lq;

              for(k=0; k<n_imm_covariate; k++)
              {
                 my_log_L_lq_imm[k] += factor1 * Q_lq_imm[k] + factor2 * Q_lq * Q_imm[k];
                 my_log_L_imm[k] += factor1 * Q_imm[k];
                 for(l=k; l<n_imm_covariate; l++)
                    my_log_L_imm_imm[k][l] += factor1 * Q_imm_imm[k][l] + factor2 * Q_imm[k] * Q_imm[l];
              }

           }


           // if people in the comunity do not share the same exposure/risk history,
           // calculate individual-level exposure/risk history.
           if(cfg_pars.common_contact_history_within_community == 0 && person->pre_immune == 0)
           {
              INITIALIZE_TEMP

              if(cfg_pars.adjust_for_left_truncation == 1 && community->earliest_idx_day_ill != MISSING)
                    start_day = max(community->day_epi_start, community->earliest_idx_day_ill - cfg_pars.max_incubation + 1);
              else  start_day = community->day_epi_start;
              day_last_obs = min(community->day_epi_stop, person->day_exit);
              day_last_obs = min(community->day_last_followup, day_last_obs);
              stop_day = (person->infection == 1)? person->day_infection_upper:day_last_obs;
              if(stop_day < start_day)
              {
                 printf("community_derivatives (line 1708): community %d, person %d, stop day %d < start day %d\n", community->id, person->id, stop_day, start_day);
                 exit(0);
              }

              for(t=start_day; t<=stop_day; t++)
              {
                 r = t - start_day;
                 rr = t - community->day_epi_start; //the time reference for risk history is community->day_epi_start, not start_day.

                 INITIALIZE_LOG_E

                 UPDATE_LOG_E_C2P
                 //The difference in infectiousness level between symptomatic and asymptomatic cases
                 //is adjusted by one additional covariate that denote the symptom status of the infective person.
                 
                 UPDATE_LOG_E_P2P

                 UPDATE_LOG_EE
                 
                 UPDATE_CUM_LOG_EE

              } /* end of day t */
           }

           if( person->pre_immune == 0)
           {
              if(cfg_pars.adjust_for_left_truncation == 1 && community->earliest_idx_day_ill != MISSING)
                    start_day = max(community->day_epi_start, community->earliest_idx_day_ill - cfg_pars.max_incubation + 1);
              else  start_day = community->day_epi_start;
              day_last_obs = min(community->day_epi_stop, person->day_exit);
              day_last_obs = min(community->day_last_followup, day_last_obs);
              stop_day = (person->infection == 1)? person->day_infection_upper:day_last_obs;

              if(stop_day < start_day)
              {
                 printf("community_derivatives (line 1743): community %d, person %d, stop day %d < start day %d\n", community->id, person->id, stop_day, start_day);
                 exit(0);
              }
              else
              {
                 INITIALIZE_L
                 INITIALIZE_LOG_L

                 /******************************************************************************************
                  calculate likelihood for a infected case with symptoms 
                  ******************************************************************************************/
                    
                 if(person->infection == 1) /* if he was infected */
                 {
                    inf1 = person->day_infection_lower;
                    inf2 = person->day_infection_upper;
                    //if(community->id == 67) printf("check 2.1: person %d  inf1=%d  inf2=%d\n", person->id, inf1, inf2);
                    if(inf2 >= start_day)
                    {
                       // it could happen that inf1 < start_day < inf2
                       inf1 = max(inf1, start_day);

                       INITIALIZE_CUM_LOG_E

                       for(t=inf1; t<=inf2; t++)
                       {
                          pdf = pdf_incubation[person->day_ill - cfg_pars.min_incubation - t];
                          r = t - start_day;
                          if(ee[r] < 0.0 || ee[r] > 1.0)
                          {
                             printf("community_derivatives (line 1772): infected person id=%d: t=%d ee[r]=%e\n", person->id, t, ee[r]);
                             error = 1;
                             exit(0);
                          }
                             
                          if(ee[r] >= 0.0 && ee[r] < 1.0)
                          {
                             UPDATE_LOG_DAY_L_FOR_INFECTION
                             // now consider pathogenicity
                             if(n_u_mode > 0)
                             {
                                rr = t - community->day_epi_start; // track covariate for pathogenicity, for which timeline starts from community->day_epi_start.
                                                                   // Covariates affecting pathogenicity are from the day of infection 
                                covariate_effect = 0.0;
                                for(k=0; k<n_pat_covariate; k++)
                                   covariate_effect += person->pat_covariate[k][rr] * coeff_pat[k];
                                U = inv_logit(lu[u_mode] + covariate_effect);
                                U_lu = U * (1 - U);
                                U_lu_lu = U * (1 - U) * (1 - 2 * U);

                                for(k=0; k<n_pat_covariate; k++)
                                {
                                   U_pat[k] = person->pat_covariate[k][rr] * U * (1 - U);
                                   U_lu_pat[k] = person->pat_covariate[k][rr] * U * (1 - U) * (1 - 2 * U);
                                   for(l=k; l<n_pat_covariate; l++)
                                      U_pat_pat[k][l] = person->pat_covariate[k][rr] * person->pat_covariate[l][rr] 
                                                      * U * (1 - U) * (1 - 2 * U);
                                }
                                pr = pdf * (person->symptom * U + (1 - person->symptom) * (1 - U));
                                pr_lu = bipow(-1, (1-person->symptom)) * pdf * U_lu;
                                pr_lu_lu = bipow(-1, (1-person->symptom)) * pdf * U_lu_lu;
                                for(k=0; k<n_pat_covariate; k++)  
                                {
                                   pr_lu_pat[k] = bipow(-1, (1-person->symptom)) * pdf * U_lu_pat[k];
                                   pr_pat[k] = bipow(-1, (1-person->symptom)) * pdf * U_pat[k];
                                   for(l=k; l<n_pat_covariate; l++)
                                      pr_pat_pat[k][l] = bipow(-1, (1-person->symptom)) * pdf * U_pat_pat[k][l];
                                }
                             }
                             else  pr = pdf;
         
                             UPDATE_L_FOR_INFECTION
                             //printf("i=%d  t=%d  L_lu=%e  L_lb_lu[0]=%e  L_lp_lu[0]=%e\n", 
                             //        person->id, t, L_lu, L_lb_lu[0], L_lp_lu[0]); 
                          }
                          UPDATE_CUM_LOG_E           
                       }
                       if(L <= 0.0)
                       {
                          printf("community_derivatives (line 1821):\n");
                          printf("Infection likelihood must be positive!  L=%f\n", L);
                          printf("Person id=%d  start_day=%d  day_ill=%d  inf1=%d  inf2=%d\n", 
                                  person->id, start_day, person->day_ill, inf1, inf2);
                          for(k=0; k<n_b_mode; k++)  printf("b[%d]=%e  ", k, b[k]);   
                          for(k=0; k<n_p_mode; k++)  printf("p[%d]=%e  ", k, p[k]);
                          printf("\n");     
                          for(t=inf1;t<=inf2;t++)
                          {
                             r = t - start_day;
                             rr = t - community->day_epi_start;
                             printf("Day %d r=%d  rr=%d:\n", t, r, rr);
                             printf("ee[]=%e 1-ee[]=%e\n", ee[r], 1.0-ee[r]);
                             
                             if(cfg_pars.common_contact_history_within_community == 1)
                             { 
                                show_c2p_risk_in_risk_class(community->id, t, t);
                                show_p2p_risk_in_risk_class(community->id, t, t);
                             }   
                             else 
                             {
                                show_c2p_risk(person, t, t);  
                                show_p2p_risk(person, t, t);
                             }   
                          } 
                          error = 1;
                          exit(0);
                       }
                       
                       UPDATE_LOG_L  
                       //if(community->id == 67) 
                       //   printf("check 2.2: i=%d  log_L_lb_lu[0]=%e  log_L_lp_lu[0]=%e\n", 
                       //           person->id, log_L_lb_lu[0], log_L_lp_lu[0]); 
                       CONCATENATE_ESCAPE_HISTORY
                       
                       //if(person->id == 3)
                       //{
                       //   show_par_effective(par_effective);
                       //   printf("day ill=%d  symptom=%d  log_L=%e  first=%e  second=%e\n", 
                       //   person->day_ill, person->symptom, log_L, log_L_lb[0], log_L_lb_lb[0][0]);
                       //   exit(0);
                       //}
                    }
                 }
                 /******************************************************************************************
                  calculate likelihood for an escaped subject, but need to consider right censoring. 
                  ******************************************************************************************/
                 else  /* if never been infected */
                 {
                       
                    if(cfg_pars.adjust_for_right_censoring == 0)
                    {
                       r = stop_day - start_day;
                       
                       log_L = cum_log_ee[r];
                       for(j=0; j<n_b_mode; j++)
                       {
                          log_L_lb[j] += cum_log_ee_lb[r][j];
                          log_L_lb_lb[j][j] += cum_log_ee_lb_lb[r][j];
                          for(k=0; k<n_c2p_covariate; k++)
                             log_L_lb_c2p[j][k] += cum_log_ee_lb_c2p[r][j][k];
                       }      

                       for(j=0; j<n_p_mode; j++)
                       {
                          log_L_lp[j] += cum_log_ee_lp[r][j];
                          log_L_lp_lp[j][j] += cum_log_ee_lp_lp[r][j];
                          for(k=0; k<n_p2p_covariate; k++)
                             log_L_lp_p2p[j][k] += cum_log_ee_lp_p2p[r][j][k];
                       }      

                       for(k=0; k<n_c2p_covariate; k++)
                       {
                          log_L_c2p[k] += cum_log_ee_c2p[r][k];
                          for(l=k; l<n_c2p_covariate; l++)
                             log_L_c2p_c2p[k][l] += cum_log_ee_c2p_c2p[r][k][l];
                       }

                       for(k=0; k<n_p2p_covariate; k++)
                       {
                          log_L_p2p[k] += cum_log_ee_p2p[r][k];
                          for(l=k; l<n_p2p_covariate; l++)
                             log_L_p2p_p2p[k][l] += cum_log_ee_p2p_p2p[r][k][l];
                       }
                      
                       //if(person->id == 2)
                       //{
                       //   show_par_effective(par_effective);
                       //   printf("day ill=%d  symptom=%d  log_L=%e  first=%e  second=%e\n", 
                       //   person->day_ill, person->symptom, log_L, log_L_lb[0], log_L_lb_lb[0][0]);
                       //   exit(0);
                       //}
                    }
                    else
                    {
                       inf1 = stop_day - cfg_pars.max_incubation + 1;
                       inf2 = stop_day - cfg_pars.min_incubation;
                       if(inf2 >= start_day)
                       {
                          // it could happen that inf1 < start_day < inf2
                          inf1 = max(inf1, start_day);

                          if(inf2 >= inf1)
                          {
                             /* potential infection during days from inf1 to inf2*/
                             INITIALIZE_CUM_LOG_E

                             for(t=inf1; t<=inf2; t++)
                             {
                                pr = sdf_incubation[inf2 - t];
                                r = t - start_day;
                                if(ee[r] < 0.0 || ee[r] > 1.0)
                                {
                                   printf("community_derivatives (line 1933): escaped person id=%d: t=%d ee[r]=%e\n", person->id, t, ee[r]);
                                   error = 1;
                                   exit(0);
                                }

                                if(ee[r] >= 0.0 && ee[r] < 1.0)
                                {
                                   UPDATE_LOG_DAY_L_FOR_INFECTION
                                   // escapes should not involve u. 
                                   // We have already set pr_lu, pr_lu_lu, pr_pat, pr_pat_pat to 0
                                   // No need to reset them here
                                   UPDATE_L_FOR_INFECTION
                                }
                                UPDATE_CUM_LOG_E
                             }

                             // escape during days from inf1 to inf2 
                             UPDATE_LOG_DAY_L_FOR_ESCAPE
                             UPDATE_L_FOR_ESCAPE

                             if(L <= 0.0)
                             {
                                printf("community_derivatives (line 1955): escaped person id=%d: L=%e is not positive\n", person->id, L);
                                error = 1;
                                exit(0);
                             }
                             UPDATE_LOG_L
                          }
                          CONCATENATE_ESCAPE_HISTORY
                          //printf("i=%d  log_L_lb_lu[0]=%e  log_L_lp_lu[0]=%e\n", 
                          //      person->id, log_L_lb_lu[0], log_L_lp_lu[0]); 
                       }
                    }                 
                 } 

                 UPDATE_MY_LOG_L
              }
           } //if(person->pre_immune == 0)

           UPDATE_SCORE_AND_INFO

           log_L_all += person->weight * my_log_L;
           ptr_integer = ptr_integer->next;
        } //while(ptr_integer != NULL)
        ptr_class = ptr_class->next;
     } //while(ptr_class != NULL)
     //if(DEBUG == 1)  printf("%d: %e %e %e\n", person->id, score_lb[0], score_lp[0], score_lu[0]);
  } /*end if(community->size > 0)*/
 

end:
  UPDATE_FIRST_AND_SECOND
  if(!(log_L_all <= 0))  log_L_all = -1e200;
  (* log_likelihood) = log_L_all;

  free(par); 
  return(error);
}


