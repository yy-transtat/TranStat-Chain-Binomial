
double community_log_likelihood_free_of_illness(COMMUNITY *community, int day_free_of_illness, double *par_effective)
{
  int i, j, h, k, l, m, n;
  int id_c2p_group, converge, error=0;
  int start_day, stop_day, max_stop_day, min_size, mean_inc;
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
  id_c2p_group = community->c2p_group;
  if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.use_index_cases_to_improve_b == 1 && community->size_idx > 0 && community->ignore == 0) 
  {
     // index cases are now included in the risk classes.
     ptr_class = community->risk_class;
     while(ptr_class != NULL)
     {
        if(cfg_pars.common_contact_history_within_community == 1)
        {
           //printf("person %d is chosen\n", person->id);
           temp = 0.0;
           // start likelihood calculation from the first infection day of the c2p phase. As we do not observe infeciton day,
           // we use earliest_idx_day_ill - mean_incubation as a proxy.
           // Note: the community risk hisotry should start from some day before this likelihood starting day.
           // Set the community->day_epi_start to a day at least as early as the phase-specific earliest_idx_day_ill - max_incubation
           start_day = community->day_epi_start;
           stop_day = community->day_epi_stop;
           //printf("start_day=%d  stop_day=%d  day_epi_stop=%d\n", start_day, stop_day, community->day_epi_stop);
           
           if(stop_day < start_day)
           {
              printf("community_log_likelihood_free_of_illness (A): community %d, stop day %d < start day %d\n", community->id, stop_day, start_day);
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
              /*ptr_risk = ptr_class->risk_history[rr].p2p_risk;
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
                    
                    if(f <= 0.0 || f > 1.0)
                    {
                       printf("Full model: f=%e\n", f);
                       error = 1;
                       goto end;
                    }
                    e *= ipow(f, ptr_risk->size);
                 }
                 ptr_risk = ptr_risk->next;
              }*/

              //update log_ee_* variables, which are related the cumulative escape probabilities
              ee[r] = e;
              temp += log(ee[r]);
              cum_log_ee[r] = temp;
           } /* end of day t */
        }
        
        ptr_integer = ptr_class->member;
        while(ptr_integer != NULL)
        {
           // unlike in the non-adjusted likelihood, index case should also contribute
           person = people + ptr_integer->id;
           //if(n_u_mode > 0)  u_mode = person->u_mode;
           //if(n_q_mode > 0)  q_mode = person->q_mode;

           my_log_L = 0.0;
           // consider prior immunity levels
           /*if(n_q_mode > 0)
           {
              covariate_effect = 0.0;
              for(k=0; k<n_imm_covariate; k++)
                 covariate_effect += person->imm_covariate[k] * coeff_imm[k];
              Q = inv_logit(lq[q_mode] + covariate_effect);

              if( person->pre_immune == 1)
                    my_log_L += log(Q);
              else
                 my_log_L += log(1 - Q);
           }*/

           //show_p2p_risk(person, 1, 5);

           // if people in the comunity do not share the same exposure/risk history,
           // calculate individual-level exposure/risk history.
           if(cfg_pars.common_contact_history_within_community == 0 && person->pre_immune == 0)
           {
              temp = 0.0;
              start_day = community->day_epi_start;
              stop_day = community->day_epi_stop;
              if(stop_day < start_day)
              {
                 printf("community_log_likelihood_free_of_illness (B): community %d, stop day %d < start day %d\n", community->id, stop_day, start_day);
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

                 /*ptr_risk = person->risk_history[rr].p2p_risk;
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
                 }    // if(person->p2p_contact_history[r].size > 0) */
                 ee[r] = e;

                 temp += log(ee[r]);
                 cum_log_ee[r] = temp;
              } /* end of day t */
           }
           // we should not include the earliest index case in the whole c2p phase
           if( person->pre_immune == 0 && (person->infection ==0 || person->day_ill > cfg_pars.c2p_group[id_c2p_group].earliest_idx_day_ill))
           {
              //start_day = max(community->day_epi_start, cfg_pars.c2p_group[id_c2p_group].earliest_idx_day_ill - cfg_pars.max_incubation + 1);
              start_day = community->day_epi_start;
              stop_day = day_free_of_illness;
              //printf("%d: start day=%d  stop day=%d\n", person->id, start_day, stop_day);
              if(stop_day < start_day)
              {
                 printf("community_log_likelihood_free_of_illness (C): community %d, stop day %d < start day %d\n", community->id, stop_day, start_day);
                 printf("c2p group earliest index case onset day= %d, theoretical likelihood start day=%d, community epidemic start day=%d\n", 
                         cfg_pars.c2p_group[id_c2p_group].earliest_idx_day_ill,
                         cfg_pars.c2p_group[id_c2p_group].earliest_idx_day_ill - cfg_pars.max_incubation + 1, community->day_epi_start);
                 exit(0);
              }   
              else
              {
                 L = log_L = 0.0;

                 /******************************************************************************************
                  calculate likelihood for an escaped subject, but need to consider right censoring. 
                  ******************************************************************************************/
                 inf1 = stop_day - cfg_pars.max_incubation + 1;
                 inf2 = stop_day - cfg_pars.min_incubation;
                 if(inf2 >= start_day)
                 {
                    // it could happen that inf1 < start_day < inf2
                    inf1 = max(inf1, start_day);
                    // When inf2>=inf1, we consider both infection plus symptom delay and escape.
                    // in2<inf1 (actually inf2=inf1-1) is also possible when min_incubation=max_incubation, 
                    // and in that case, we are sure escape occured up to inf1-1,
                    // and we only need the CONCATENATE_ESCAPE_HISTORY part. For more details,
                    // see the comments in function community_derivatives_IdxAdjust_denominator(). 
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
                             printf("community_loglikelihood_free_of_illness, escaped person id=%d: t=%d  ee[r]=%e\n", person->id, t, ee[r]);
                             error = 1;
                             goto end;
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
                          printf("community_loglikelihood_free_of_illness, escaped person id=%d: L=%e is not positive\n", person->id, L);
                          error = 1;
                          goto end;
                       }

                       log_L = log(L);
                    }
                    r = inf1 - start_day;
                    if(r > 0) log_L += cum_log_ee[r-1];
                 }
                 my_log_L += log_L;
              }
           } //if( person->pre_immune == 0)
           log_L_all += person->weight * my_log_L;
           //if(community->id == 42)  printf("person %d: my_log_L=%e  log_L_all=%e\n", person->id, my_log_L, log_L_all);
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


double community_log_likelihood_IdxAdjust_numerator(COMMUNITY *community, double *par_effective)
{
   double log_L_all, log_L_1, log_L_2;
   log_L_1 = community_log_likelihood_free_of_illness(community, community->earliest_idx_day_ill, par_effective);
   log_L_2 = community_log_likelihood_free_of_illness(community, community->earliest_idx_day_ill-1, par_effective);
   log_L_all = log(exp(log_L_2) - exp(log_L_1));
   if(log_L_1 >= log_L_2)
   {
      printf("community_log_likelihood_IdxAdjust_numerator: log_L_1(%e) >= log_L_2(%e)\n", log_L_1, log_L_2);
      exit(0);
   }   
   return(log_L_all);
}


double community_log_likelihood_IdxAdjust_denominator(COMMUNITY *community, double *par_effective)
{
  double log_L_all, log_L_1, log_L_2;
  int id_c2p_group;
  id_c2p_group = community->c2p_group;
  log_L_1 = community_log_likelihood_free_of_illness(community, community->day_epi_stop, par_effective);
  log_L_2 = community_log_likelihood_free_of_illness(community, cfg_pars.c2p_group[id_c2p_group].earliest_idx_day_ill, par_effective);
  log_L_all = log(exp(log_L_2) - exp(log_L_1));
   if(log_L_1 >= log_L_2)
   {
      printf("community_log_likelihood_IdxAdjust_denominator: log_L_1(%e) >= log_L_2(%e)\n", log_L_1, log_L_2);
      exit(0);
   }   
  return(log_L_all);
}




int community_derivatives_free_of_illness(COMMUNITY *community, int day_free_of_illness,  double *par_effective, double *log_likelihood, MATRIX *first, MATRIX *second)
{
  int i, j, h, k, l, m, n;
  int id_c2p_group, converge, error=0;
  int start_day, stop_day, max_stop_day, min_size;
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

  id_c2p_group = community->c2p_group;
  if(cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.use_index_cases_to_improve_b == 1 && community->size_idx > 0 && community->ignore == 0) 
  {

     ptr_class = community->risk_class;
     while(ptr_class != NULL)
     {
        if(cfg_pars.common_contact_history_within_community == 1)
        {
           INITIALIZE_TEMP

           start_day = community->day_epi_start;
           stop_day = community->day_epi_stop;
           //printf("%d: start_day=%d  stop_day=%d\n", h, start_day, stop_day);
           if(stop_day < start_day)
           {
              printf("community_derivatives_free_of_illness (A): community %d, stop day %d < start day %d\n", community->id, stop_day, start_day);
              error = 1;
              goto end;
           }   
           for(t=start_day; t<=stop_day; t++)
           {
              r = t - start_day;
              rr = t - community->day_epi_start; //the time reference for risk history is community->day_epi_start, not start_day.
                                                 //
              INITIALIZE_LOG_E

              UPDATE_LOG_E_C2P
              //The difference in infectiousness level between symptomatic and asymptomatic cases
              //is adjusted by one additional covariate that denote the symptom status of the infective person.
              
              //UPDATE_LOG_E_P2P
              
              UPDATE_LOG_EE
              
              UPDATE_CUM_LOG_EE
              //printf("t=%d  log_ee_lb=%e  log_ee_lb_lb=%e  log_ee_lp=%e  log_ee_lp_lp=%e\n", t, log_ee_lb[r][0], log_ee_lb_lb[r][0], log_ee_lp[r][0], log_ee_lp_lp[r][0]);
           } /* end of day t */
        }
        
        ptr_integer = ptr_class->member;
        while(ptr_integer != NULL)
        {
           // unlike in the non-adjusted likelihood, index case should also contribute
           person = people + ptr_integer->id;
     
           if(n_u_mode > 0)  u_mode = person->u_mode;
           if(n_q_mode > 0)  q_mode = person->q_mode;

           INITIALIZE_MY_LOG_L


           // if people in the comunity do not share the same exposure/risk history,
           // calculate individual-level exposure/risk history.
           if(cfg_pars.common_contact_history_within_community == 0 && person->pre_immune == 0)
           {
              INITIALIZE_TEMP

              start_day = community->day_epi_start;
              stop_day = community->day_epi_stop;
              if(stop_day < start_day)
              {
                 printf("community_derivatives_free_of_illness (B): community %d, stop day %d < start day %d\n", community->id, stop_day, start_day);
                 error = 1;
                 goto end;
              }   
              for(t=start_day; t<=stop_day; t++)
              {
                 r = t - start_day;
                 rr = t - community->day_epi_start; //the time reference for risk history is community->day_epi_start, not start_day.

                 INITIALIZE_LOG_E

                 UPDATE_LOG_E_C2P
                 //The difference in infectiousness level between symptomatic and asymptomatic cases
                 //is adjusted by one additional covariate that denote the symptom status of the infective person.
                 
                 //UPDATE_LOG_E_P2P
                 
                 UPDATE_LOG_EE
                 
                 UPDATE_CUM_LOG_EE
              } /* end of day t */
           }

           if( person->pre_immune == 0 && (person->infection ==0 || person->day_ill > cfg_pars.c2p_group[id_c2p_group].earliest_idx_day_ill))
           {
              //start_day = max(community->day_epi_start, cfg_pars.c2p_group[id_c2p_group].earliest_idx_day_ill - cfg_pars.max_incubation + 1);
              start_day = community->day_epi_start;
              stop_day = day_free_of_illness;

              if(stop_day < start_day)
              {
                 printf("community_derivatives_free_of_illness (C): community %d, stop day %d < start day %d\n", community->id, stop_day, start_day);
                 printf("c2p group earliest index case onset day= %d, theoretical likelihood start day=%d, community epidemic start day=%d\n", 
                         cfg_pars.c2p_group[id_c2p_group].earliest_idx_day_ill,
                         cfg_pars.c2p_group[id_c2p_group].earliest_idx_day_ill - cfg_pars.max_incubation + 1, community->day_epi_start);
                 error = 1;
                 goto end;
              }   
              else
              {
                 INITIALIZE_L
                 INITIALIZE_LOG_L

                 /******************************************************************************************
                  calculate likelihood for an escaped subject, but need to consider right censoring. 
                  ******************************************************************************************/
                 inf1 = stop_day - cfg_pars.max_incubation + 1;
                 inf2 = stop_day - cfg_pars.min_incubation;
                 if(inf2 >= start_day)
                 {
                    // it could happen that inf1 < start_day < inf2
                    inf1 = max(inf1, start_day);

                    // When inf2>=inf1, we consider both infection plus symptom delay and escape.
                    // in2<inf1 (actually inf2=inf1-1) is also possible when min_incubation=max_incubation, 
                    // and in that case, we are sure escape occured up to inf1-1,
                    // and we only need the CONCATENATE_ESCAPE_HISTORY part. For more details,
                    // see the comments in function community_derivatives_IdxAdjust_denominator(). 
                    if(inf2 >= inf1)
                    {
                       INITIALIZE_CUM_LOG_E

                       // potential infection during days from inf1 to inf2
                       for(t=inf1; t<=inf2; t++)
                       {
                          pr = sdf_incubation[inf2 - t];
                          r = t - start_day;
                          if(ee[r] < 0.0 || ee[r] > 1.0)
                          {
                             printf("community_derivatives_free_of_illness, escaped person id=%d: t=%d ee[r]=%e\n", person->id, t, ee[r]);
                             error = 1;
                             goto end;
                          }
 
                          if(ee[r] >= 0.0 && ee[r] < 1.0)
                          {
                             UPDATE_LOG_DAY_L_FOR_INFECTION
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
                          printf("community_derivatives_free_of_illness, escaped person id=%d: L=%e is not positive\n", person->id, L);
                          error = 1;
                          goto end;
                       }

                       UPDATE_LOG_L
                    }
                    CONCATENATE_ESCAPE_HISTORY
                    //printf("i=%d  log_L_lb_lu[0]=%e  log_L_lp_lu[0]=%e\n", person->id, log_L_lb_lu[0], log_L_lp_lu[0]); 
                 }

                 UPDATE_MY_LOG_L
              }
           } //if(person->pre_immune == 0)

           UPDATE_SCORE_AND_INFO

           log_L_all += person->weight * my_log_L;
           //if(community->id == 42)  printf("person %d: my_log_L=%e  log_L_all=%e\n", person->id, my_log_L, log_L_all);
           ptr_integer = ptr_integer->next;
           //printf("Numerator person %d: %e %e %e\n", person->id, my_log_L, my_log_L_lb[0], my_log_L_lp[0]);
        } //while(ptr_integer != NULL)
        ptr_class = ptr_class->next;
     } //while(ptr_class != NULL)
  } /*end if(community->size > 0)*/
 
  UPDATE_FIRST_AND_SECOND 

end:
  if(!(log_L_all <= 0)) error = 1;
  if(error == 1)  log_L_all = -1e200;
  (* log_likelihood) = log_L_all;

  free(par); 
  return(error);
}


int community_derivatives_IdxAdjust_numerator(COMMUNITY *community, double *par_effective, double *log_likelihood, MATRIX *first, MATRIX *second)
{
   int i, j, n_par, error, error1, error2;
   double L, L_1, L_2, log_L_1, log_L_2;
   MATRIX d_L, dd_L, first_1, first_2, second_1, second_2;
 
   initialize_matrix(&d_L);
   initialize_matrix(&dd_L);
   initialize_matrix(&first_1);
   initialize_matrix(&first_2);
   initialize_matrix(&second_1);
   initialize_matrix(&second_2);
 
   n_par = cfg_pars.n_par;
   inflate_matrix(&d_L, n_par, 1, 0);
   inflate_matrix(&dd_L, n_par, n_par, 0);
   inflate_matrix(&first_1, n_par, 1, 0);
   inflate_matrix(&first_2, n_par, 1, 0);
   inflate_matrix(&second_1, n_par, n_par, 0);
   inflate_matrix(&second_2, n_par, n_par, 0);

   //if(community->id == 162)  printf("first call: day_free_of_illness=%d\n", community->earliest_idx_day_ill);
   error1 = community_derivatives_free_of_illness(community, community->earliest_idx_day_ill,  par_effective, &log_L_1, &first_1, &second_1);
   //if(community->id == 162)  printf("second call: day_free_of_illness=%d\n", community->earliest_idx_day_ill - 1);
   error2 = community_derivatives_free_of_illness(community, community->earliest_idx_day_ill-1,  par_effective, &log_L_2, &first_2, &second_2);
   error = (error1 == 1 || error2 == 1)? 1:0;

   if(error1 == 1)
   {
      log_L_1 = -1e200;
      reset(&first_1, 0);
      reset(&second_1, 0);
   }   
   if(error2 == 1)
   {
      log_L_2 = -1e200;
      reset(&first_2, 0);
      reset(&second_2, 0);
   }   
   if(log_L_1 >= log_L_2)
   {
      error = 1;
      printf("community_derivatives_IdxAdjust_numerator: log_L_1(%e) >= log_L_2(%e)\n", log_L_1, log_L_2);
      exit(0);
   }   

   L_2 = exp(log_L_2);
   L_1 = exp(log_L_1);
   L = L_2 - L_1;
   
   // d_L is the first derivatives of L; dd_L is the second serivatives of L
   (* log_likelihood) = log(L);
   for(i=0; i<n_par; i++)
   {
      d_L.data[i][0] = L_2 * first_2.data[i][0] - L_1 * first_1.data[i][0];
      first->data[i][0] = d_L.data[i][0] / L;
      for(j=0; j<=i; j++)
      {
         dd_L.data[i][j] = L_2 * (first_2.data[i][0] * first_2.data[j][0] + second_2.data[i][j]) -
                           L_1 * (first_1.data[i][0] * first_1.data[j][0] + second_1.data[i][j]);
         second->data[i][j] = - d_L.data[i][0] * d_L.data[j][0] / (L*L) + dd_L.data[i][j] / L;
         second->data[j][i] = second->data[i][j];
      }    
   }

   deflate_matrix(&d_L);
   deflate_matrix(&dd_L);
   deflate_matrix(&first_1);
   deflate_matrix(&first_2);
   deflate_matrix(&second_1);
   deflate_matrix(&second_2);
   return(error);
}





int community_derivatives_IdxAdjust_denominator(COMMUNITY *community, double *par_effective, double *log_likelihood, MATRIX *first, MATRIX *second)
{
   int i, j, n_par, id_c2p_group, error, error1, error2; 
   double L, L_1, L_2, log_L_1, log_L_2;
   MATRIX d_L, dd_L, first_1, first_2, second_1, second_2;

   n_par = cfg_pars.n_par;

   initialize_matrix(&d_L);
   initialize_matrix(&dd_L);
   initialize_matrix(&first_1);
   initialize_matrix(&first_2);
   initialize_matrix(&second_1);
   initialize_matrix(&second_2);
 
   inflate_matrix(&d_L, n_par, 1, 0);
   inflate_matrix(&dd_L, n_par, n_par, 0);
   inflate_matrix(&first_1, n_par, 1, 0);
   inflate_matrix(&first_2, n_par, 1, 0);
   inflate_matrix(&second_1, n_par, n_par, 0);
   inflate_matrix(&second_2, n_par, n_par, 0);

   id_c2p_group = community->c2p_group;
   error1 = community_derivatives_free_of_illness(community, community->day_epi_stop, par_effective, &log_L_1, &first_1, &second_1);
   error2 = community_derivatives_free_of_illness(community, cfg_pars.c2p_group[id_c2p_group].earliest_idx_day_ill,  par_effective, &log_L_2, &first_2, &second_2);
   error = (error1 == 1 || error2 == 1)? 1:0;

   if(error1 == 1)
   {
      log_L_1 = -1e200;
      reset(&first_1, 0);
      reset(&second_1, 0);
   }   
   if(error2 == 1)
   {
      log_L_2 = -1e200;
      reset(&first_2, 0);
      reset(&second_2, 0);
   }   
   if(log_L_1 >= log_L_2)
   {
      error = 1;
      printf("community_derivatives_IdxAdjust_numerator: log_L_1(%e) >= log_L_2(%e)\n", log_L_1, log_L_2);
      exit(0);
   }   
   /*if(community->id == 2)
   {
      printf("day_epi_start=%d  day_epie_stop=%d  c2p earliest idx day ill=%d\n", community->day_epi_start, community->day_epi_stop, cfg_pars.c2p_group[id_c2p_group].earliest_idx_day_ill);
      printf("error1=%d  error2=%d\n", error1, error2);
      printf("log_L_1=%e  log_L_2=%e\n", log_L_1, log_L_2);
      printf("first 1\n");
      fmprintf(&first_1);
      printf("second 1\n");
      fmprintf(&second_1);
      printf("first 2\n");
      fmprintf(&first_2);
      printf("second 2\n");
      fmprintf(&second_2);
   } */  
   L_2 = exp(log_L_2);
   L_1 = exp(log_L_1);

   L = L_2 - L_1;
   
   // d_L is the first derivatives of L; dd_L is the second serivatives of L
   (* log_likelihood) = log(L);
   for(i=0; i<n_par; i++)
   {
      d_L.data[i][0] = L_2 * first_2.data[i][0] - L_1 * first_1.data[i][0];
      first->data[i][0] = d_L.data[i][0] / L;
      for(j=0; j<=i; j++)
      {
         dd_L.data[i][j] = L_2 * (first_2.data[i][0] * first_2.data[j][0] + second_2.data[i][j]) -
                           L_1 * (first_1.data[i][0] * first_1.data[j][0] + second_1.data[i][j]);
         second->data[i][j] = - d_L.data[i][0] * d_L.data[j][0] / (L*L) + dd_L.data[i][j] / L;
         second->data[j][i] = second->data[i][j];
      }    
   }

   deflate_matrix(&d_L);
   deflate_matrix(&dd_L);
   deflate_matrix(&first_1);
   deflate_matrix(&first_2);
   deflate_matrix(&second_1);
   deflate_matrix(&second_2);
   return(error);
}


