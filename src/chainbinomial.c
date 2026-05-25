/* chainbinomial.c
 * Unified compilation unit for the ChainBinomial R package.
 * Includes all shared headers exactly once, declares global state,
 * and provides R-callable wrappers for simulation and estimation.
 *
 * Author: Yang Yang
 * Department of Biostatistics and Emerging Pathogens Institute
 * University of Florida
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

/* Macros (matches main.c exactly; duplicates removed) */
#define logit(x)     ( log( (x) / ( 1.0 - (x) ) ) )
#define inv_logit(x) ( 1.0 / ( 1.0 + exp(-(x)) ) )
#define bipow(x,y)   (((y)==0)? 1:(x))
#define max(x, y)    ((x)>(y)?(x):(y))
#define min(x, y)    ((x)<(y)?(x):(y))
#define square(x)    ( (x) * (x) )
#define MISSING             99999
#define INFINITY_INTEGER    1000000
#define close_to_0          1e-20
#define close_to_1          0.99999
#define close_to_05         0.49999
#define num_ini_estimation  50
#define num_ini_test        10

/* Shared library headers (included once here; not in main.c or pop_gen_household.c) */
#include "matrix.h"
#include "mathfunc.h"
#include "distribution.h"
#include "statistics.h"
#include "optimization.h"

/* Data structures */
#include "datastruct.h"

/* Global state required by the estimation headers below */
CFG_PARS  cfg_pars;
PEOPLE   *people;
int       p_size;
COMMUNITY *community;
int       n_community;
int       DEBUG = 0;
int       max_epi_duration;
long      seed;
int       serial_number;
double   *pdf_incubation, *sdf_incubation;
int       size_sample_states;
double   *importance_weight;
int       ITER;
int       n_need_MCEM, n_need_OEM;
long      TIME_SAMPLING_SEC, TIME_SAMPLING_CPU;
long      TIME_ITERATION_SEC, TIME_ITERATION_CPU;
long      TIME_VARIANCE_SEC, TIME_VARIANCE_CPU;

/* Forward declaration required because core.h calls get_size() before it is defined.
   The original standalone build relied on C89 implicit declarations; strict C99 requires this. */
int get_size(FILE *in);

/* Estimation and simulation logic (headers use the globals declared above) */
#include "config.h"
#include "utility.h"
#include "simulate.h"
#include "preamble.h"
#include "derivatives.h"
#include "derivatives_IdxAdjust.h"
#include "estimation.h"
#include "test_cs_hh_improve.h"
#include "core.h"

int get_size(FILE *in)
{
    int c, i = 0;
    rewind(in);
    while ((c = getc(in)) != EOF)
        if (c == '\n') i++;
    return i;
}

/* =========================================================
 * R interface
 * ========================================================= */

/* Several of our macros conflict with identifiers in R's Rinternals.h.
   Undefine them here; all C code above has already consumed them. */
#undef max
#undef min
#undef MISSING          /* R's Rinternals.h declares int (MISSING)(SEXP x) */
#define CB_MISSING 99999 /* replacement sentinel for use in the R interface below */
#undef INFINITY_INTEGER

/* Replacements for the max/min macros that were undefined above.
   Used only within the R wrapper functions below. */
#define CB_max(a, b) ((a) > (b) ? (a) : (b))
#define CB_min(a, b) ((a) < (b) ? (a) : (b))

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* ---------------------------------------------------------
 * r_gen_population
 *
 * Generate a pseudo-population for chain binomial simulation.
 * Mirrors the logic of pop_gen_household.c::main() but returns
 * data directly to R as a named list of three data frames.
 *
 * Arguments (all SEXP wrappers of scalars):
 *   r_n_community          – number of communities
 *   r_community_size       – individuals per community
 *   r_day_epi_stop         – last day of epidemic
 *   r_case_ascertained     – 1 = first member is pre-seeded index case
 *   r_cluster_randomization– 1 = cluster-level treatment, 0 = individual
 *   r_seed                 – RNG seed
 *   r_prop_idx             – proportion of index cases untreated
 *   r_prop_contact         – proportion of contacts untreated
 *
 * Returns a list: $pop, $community, $time_dep_covariate (data frames)
 * --------------------------------------------------------- */
SEXP r_gen_population(SEXP r_n_community, SEXP r_community_size,
                      SEXP r_day_epi_stop, SEXP r_case_ascertained,
                      SEXP r_cluster_randomization, SEXP r_seed,
                      SEXP r_prop_idx, SEXP r_prop_contact)
{
    int h, i, j, k, l, t, col, row;
    int loc_p_size, loc_n_community, community_size;
    int day_epi_start, day_epi_stop;
    int case_ascertained, cluster_randomization, n_epi_days, n_tdcov_rows;
    double prop_idx, prop_contact;
    PEOPLE    *loc_people, *person;
    COMMUNITY *loc_community;

    /* Unpack R arguments */
    loc_n_community      = INTEGER(r_n_community)[0];
    community_size       = INTEGER(r_community_size)[0];
    day_epi_start        = 1;
    day_epi_stop         = INTEGER(r_day_epi_stop)[0];
    case_ascertained     = INTEGER(r_case_ascertained)[0];
    cluster_randomization= INTEGER(r_cluster_randomization)[0];
    prop_idx             = REAL(r_prop_idx)[0];
    prop_contact         = REAL(r_prop_contact)[0];

    seed = (long)INTEGER(r_seed)[0];
    mt_init(seed);

    loc_p_size = loc_n_community * community_size;

    /* Allocate population */
    loc_people = (PEOPLE *)malloc((size_t)(loc_p_size * sizeof(PEOPLE)));
    if (!loc_people) error("Cannot allocate memory for people array");

    for (i = 0, person = loc_people; i < loc_p_size; i++, person++) {
        person->id        = i;
        person->community = i / community_size;
        make_1d_array_double(&person->time_ind_covariate, 1, 0.0);
        make_2d_array_double(&person->time_dep_covariate, 1,
                             day_epi_stop - day_epi_start + 1, 0.0);
        person->pre_immune = 0;
        person->infection  = 0;
        person->symptom    = 0;
        person->day_ill    = CB_MISSING;
        person->exit       = 0;
        person->day_exit   = CB_MISSING;
        person->idx        = 0;
        person->u_mode     = 0;
        person->q_mode     = 0;
        person->ignore     = 0;
        person->weight     = 1;
    }

    /* Allocate community structure */
    loc_community = (COMMUNITY *)malloc((size_t)(loc_n_community * sizeof(COMMUNITY)));
    if (!loc_community) { free(loc_people); error("Cannot allocate memory for community array"); }

    for (h = 0; h < loc_n_community; h++) {
        loc_community[h].size              = 0;
        loc_community[h].member            = (int *)malloc((size_t)(community_size * sizeof(int)));
        loc_community[h].day_epi_start     = day_epi_start;
        loc_community[h].day_epi_stop      = day_epi_stop;
        loc_community[h].day_last_followup = day_epi_stop;
        loc_community[h].c2p_group         = 0;
    }

    for (i = 0, person = loc_people; i < loc_p_size; i++, person++) {
        h = person->community;
        loc_community[h].member[loc_community[h].size] = i;
        loc_community[h].size++;
    }

    /* Seed index cases in case-ascertained design */
    if (case_ascertained == 1) {
        for (h = 0; h < loc_n_community; h++) {
            i = loc_community[h].member[0];
            loc_people[i].idx       = 1;
            loc_people[i].infection = 1;
            loc_people[i].symptom   = 1;
            loc_people[i].day_ill   = 1;
        }
    }

    /* Assign time-independent treatment covariate.
       prop_idx    = P(covariate=1 | index case)
       prop_contact = P(covariate=1 | contact) */
    if (cluster_randomization == 1) {
        for (h = 0; h < loc_n_community; h++) {
            k = (runiform(&seed) < prop_idx);
            l = (runiform(&seed) < prop_contact);
            for (j = 0; j < loc_community[h].size; j++) {
                i = loc_community[h].member[j];
                if (loc_people[i].idx == 1) loc_people[i].time_ind_covariate[0] = k;
                else                        loc_people[i].time_ind_covariate[0] = l;
            }
        }
    } else {
        for (h = 0; h < loc_n_community; h++) {
            for (j = 0; j < loc_community[h].size; j++) {
                i = loc_community[h].member[j];
                if (loc_people[i].idx == 1)
                    loc_people[i].time_ind_covariate[0] = (runiform(&seed) < prop_idx);
                else
                    loc_people[i].time_ind_covariate[0] = (runiform(&seed) < prop_contact);
            }
        }
    }

    /* Generate time-dependent covariate from standard normal */
    for (i = 0; i < loc_p_size; i++)
        for (t = day_epi_start; t <= day_epi_stop; t++)
            loc_people[i].time_dep_covariate[0][t-1] = rnorm(0.0, 1.0, &seed);

    /* ---- Build R return value ---- */
    n_epi_days    = day_epi_stop - day_epi_start + 1;
    n_tdcov_rows  = loc_p_size * n_epi_days;

    /* pop vectors: community/disease info only (13 columns, no covariate) */
    SEXP pop_id, pop_comm, pop_pre_immune, pop_infection, pop_symptom;
    SEXP pop_day_ill, pop_exit, pop_day_exit, pop_idx, pop_u_mode, pop_q_mode;
    SEXP pop_weight, pop_ignore;
    /* time-independent covariate: separate data frame (id, value) */
    SEXP tic_id, tic_val;
    /* community vectors */
    SEXP com_id, com_ds, com_dp, com_dfu, com_c2p;
    /* time-dependent covariate vectors */
    SEXP tdc_id, tdc_day_start, tdc_day_stop, tdc_val;
    /* data frame and result objects */
    SEXP pop_df, pop_nm, tic_df, tic_nm, com_df, com_nm, tdc_df, tdc_nm;
    SEXP result, res_nm, rn;

    /* pop vectors (1-13) */
    PROTECT(pop_id        = allocVector(INTSXP,  loc_p_size));
    PROTECT(pop_comm      = allocVector(INTSXP,  loc_p_size));
    PROTECT(pop_pre_immune= allocVector(INTSXP,  loc_p_size));
    PROTECT(pop_infection = allocVector(INTSXP,  loc_p_size));
    PROTECT(pop_symptom   = allocVector(INTSXP,  loc_p_size));
    PROTECT(pop_day_ill   = allocVector(INTSXP,  loc_p_size));
    PROTECT(pop_exit      = allocVector(INTSXP,  loc_p_size));
    PROTECT(pop_day_exit  = allocVector(INTSXP,  loc_p_size));
    PROTECT(pop_idx       = allocVector(INTSXP,  loc_p_size));
    PROTECT(pop_u_mode    = allocVector(INTSXP,  loc_p_size));
    PROTECT(pop_q_mode    = allocVector(INTSXP,  loc_p_size));
    PROTECT(pop_weight    = allocVector(REALSXP, loc_p_size));
    PROTECT(pop_ignore    = allocVector(INTSXP,  loc_p_size)); /* 13 */

    /* time-independent covariate vectors (14-15) */
    PROTECT(tic_id  = allocVector(INTSXP,  loc_p_size));
    PROTECT(tic_val = allocVector(REALSXP, loc_p_size)); /* 15 */

    /* community vectors (16-20) */
    PROTECT(com_id  = allocVector(INTSXP, loc_n_community));
    PROTECT(com_ds  = allocVector(INTSXP, loc_n_community));
    PROTECT(com_dp  = allocVector(INTSXP, loc_n_community));
    PROTECT(com_dfu = allocVector(INTSXP, loc_n_community));
    PROTECT(com_c2p = allocVector(INTSXP, loc_n_community)); /* 20 */

    /* time-dependent covariate vectors (21-24) */
    PROTECT(tdc_id        = allocVector(INTSXP,  n_tdcov_rows));
    PROTECT(tdc_day_start = allocVector(INTSXP,  n_tdcov_rows));
    PROTECT(tdc_day_stop  = allocVector(INTSXP,  n_tdcov_rows));
    PROTECT(tdc_val       = allocVector(REALSXP, n_tdcov_rows)); /* 24 */

    /* Fill pop (disease/community info only) */
    for (i = 0; i < loc_p_size; i++) {
        INTEGER(pop_id)[i]         = loc_people[i].id;
        INTEGER(pop_comm)[i]       = loc_people[i].community;
        INTEGER(pop_pre_immune)[i] = loc_people[i].pre_immune;
        INTEGER(pop_infection)[i]  = loc_people[i].infection;
        INTEGER(pop_symptom)[i]    = loc_people[i].symptom;
        INTEGER(pop_day_ill)[i]    = loc_people[i].day_ill;
        INTEGER(pop_exit)[i]       = loc_people[i].exit;
        INTEGER(pop_day_exit)[i]   = loc_people[i].day_exit;
        INTEGER(pop_idx)[i]        = loc_people[i].idx;
        INTEGER(pop_u_mode)[i]     = loc_people[i].u_mode;
        INTEGER(pop_q_mode)[i]     = loc_people[i].q_mode;
        REAL(pop_weight)[i]        = loc_people[i].weight;
        INTEGER(pop_ignore)[i]     = loc_people[i].ignore;
    }

    /* Fill time-independent covariate */
    for (i = 0; i < loc_p_size; i++) {
        INTEGER(tic_id)[i] = loc_people[i].id;
        REAL(tic_val)[i]   = loc_people[i].time_ind_covariate[0];
    }

    /* Fill community */
    for (h = 0; h < loc_n_community; h++) {
        INTEGER(com_id)[h]  = h;
        INTEGER(com_ds)[h]  = loc_community[h].day_epi_start;
        INTEGER(com_dp)[h]  = loc_community[h].day_epi_stop;
        INTEGER(com_dfu)[h] = loc_community[h].day_last_followup;
        INTEGER(com_c2p)[h] = loc_community[h].c2p_group;
    }

    /* Fill time-dependent covariate (interval format: id, day_start, day_stop, value).
       For this simulation each day is drawn independently so day_start == day_stop == t.
       The interval format generalises to covariates that are constant across segments
       (e.g., vaccination status: 0 before day T, 1 on and after day T). */
    row = 0;
    for (i = 0; i < loc_p_size; i++) {
        for (t = day_epi_start; t <= day_epi_stop; t++) {
            INTEGER(tdc_id)[row]        = i;
            INTEGER(tdc_day_start)[row] = t;
            INTEGER(tdc_day_stop)[row]  = t;
            REAL(tdc_val)[row]          = loc_people[i].time_dep_covariate[0][t-1];
            row++;
        }
    }

    /* Free C memory before building R objects */
    for (i = 0; i < loc_p_size; i++) {
        free(loc_people[i].time_ind_covariate);
        free_2d_array_double(loc_people[i].time_dep_covariate);
    }
    free(loc_people);
    for (h = 0; h < loc_n_community; h++) free(loc_community[h].member);
    free(loc_community);

    /* Build pop data frame — 13 columns, community/disease info only (24, 25) */
    PROTECT(pop_df = allocVector(VECSXP, 13)); /* 24 */
    PROTECT(pop_nm = allocVector(STRSXP, 13)); /* 25 */
    SET_VECTOR_ELT(pop_df,  0, pop_id);        SET_STRING_ELT(pop_nm,  0, mkChar("id"));
    SET_VECTOR_ELT(pop_df,  1, pop_comm);      SET_STRING_ELT(pop_nm,  1, mkChar("community"));
    SET_VECTOR_ELT(pop_df,  2, pop_pre_immune);SET_STRING_ELT(pop_nm,  2, mkChar("pre_immune"));
    SET_VECTOR_ELT(pop_df,  3, pop_infection); SET_STRING_ELT(pop_nm,  3, mkChar("infection"));
    SET_VECTOR_ELT(pop_df,  4, pop_symptom);   SET_STRING_ELT(pop_nm,  4, mkChar("symptom"));
    SET_VECTOR_ELT(pop_df,  5, pop_day_ill);   SET_STRING_ELT(pop_nm,  5, mkChar("day_ill"));
    SET_VECTOR_ELT(pop_df,  6, pop_exit);      SET_STRING_ELT(pop_nm,  6, mkChar("exit"));
    SET_VECTOR_ELT(pop_df,  7, pop_day_exit);  SET_STRING_ELT(pop_nm,  7, mkChar("day_exit"));
    SET_VECTOR_ELT(pop_df,  8, pop_idx);       SET_STRING_ELT(pop_nm,  8, mkChar("idx"));
    SET_VECTOR_ELT(pop_df,  9, pop_u_mode);    SET_STRING_ELT(pop_nm,  9, mkChar("u_mode"));
    SET_VECTOR_ELT(pop_df, 10, pop_q_mode);    SET_STRING_ELT(pop_nm, 10, mkChar("q_mode"));
    SET_VECTOR_ELT(pop_df, 11, pop_weight);    SET_STRING_ELT(pop_nm, 11, mkChar("weight"));
    SET_VECTOR_ELT(pop_df, 12, pop_ignore);    SET_STRING_ELT(pop_nm, 12, mkChar("ignore"));
    setAttrib(pop_df, R_NamesSymbol, pop_nm);
    setAttrib(pop_df, R_ClassSymbol, mkString("data.frame"));
    PROTECT(rn = allocVector(INTSXP, 2)); /* 26 */
    INTEGER(rn)[0] = NA_INTEGER; INTEGER(rn)[1] = -loc_p_size;
    setAttrib(pop_df, R_RowNamesSymbol, rn);
    UNPROTECT(1); /* rn; stack = 25 */

    /* Build time_ind_covariate data frame — 2 columns: id, value (26, 27) */
    PROTECT(tic_df = allocVector(VECSXP, 2)); /* 26 */
    PROTECT(tic_nm = allocVector(STRSXP, 2)); /* 27 */
    SET_VECTOR_ELT(tic_df, 0, tic_id);  SET_STRING_ELT(tic_nm, 0, mkChar("id"));
    SET_VECTOR_ELT(tic_df, 1, tic_val); SET_STRING_ELT(tic_nm, 1, mkChar("value"));
    setAttrib(tic_df, R_NamesSymbol, tic_nm);
    setAttrib(tic_df, R_ClassSymbol, mkString("data.frame"));
    PROTECT(rn = allocVector(INTSXP, 2)); /* 28 */
    INTEGER(rn)[0] = NA_INTEGER; INTEGER(rn)[1] = -loc_p_size;
    setAttrib(tic_df, R_RowNamesSymbol, rn);
    UNPROTECT(1); /* rn; stack = 27 */

    /* Build community data frame (28, 29) */
    PROTECT(com_df = allocVector(VECSXP, 5)); /* 28 */
    PROTECT(com_nm = allocVector(STRSXP, 5)); /* 29 */
    SET_VECTOR_ELT(com_df, 0, com_id);  SET_STRING_ELT(com_nm, 0, mkChar("community"));
    SET_VECTOR_ELT(com_df, 1, com_ds);  SET_STRING_ELT(com_nm, 1, mkChar("day_epi_start"));
    SET_VECTOR_ELT(com_df, 2, com_dp);  SET_STRING_ELT(com_nm, 2, mkChar("day_epi_stop"));
    SET_VECTOR_ELT(com_df, 3, com_dfu); SET_STRING_ELT(com_nm, 3, mkChar("day_last_followup"));
    SET_VECTOR_ELT(com_df, 4, com_c2p); SET_STRING_ELT(com_nm, 4, mkChar("c2p_group"));
    setAttrib(com_df, R_NamesSymbol, com_nm);
    setAttrib(com_df, R_ClassSymbol, mkString("data.frame"));
    PROTECT(rn = allocVector(INTSXP, 2)); /* 30 */
    INTEGER(rn)[0] = NA_INTEGER; INTEGER(rn)[1] = -loc_n_community;
    setAttrib(com_df, R_RowNamesSymbol, rn);
    UNPROTECT(1); /* rn; stack = 29 */

    /* Build time_dep_covariate data frame — 4 columns (31, 32) */
    PROTECT(tdc_df = allocVector(VECSXP, 4)); /* 31 */
    PROTECT(tdc_nm = allocVector(STRSXP, 4)); /* 32 */
    SET_VECTOR_ELT(tdc_df, 0, tdc_id);        SET_STRING_ELT(tdc_nm, 0, mkChar("id"));
    SET_VECTOR_ELT(tdc_df, 1, tdc_day_start); SET_STRING_ELT(tdc_nm, 1, mkChar("day_start"));
    SET_VECTOR_ELT(tdc_df, 2, tdc_day_stop);  SET_STRING_ELT(tdc_nm, 2, mkChar("day_stop"));
    SET_VECTOR_ELT(tdc_df, 3, tdc_val);       SET_STRING_ELT(tdc_nm, 3, mkChar("value"));
    setAttrib(tdc_df, R_NamesSymbol, tdc_nm);
    setAttrib(tdc_df, R_ClassSymbol, mkString("data.frame"));
    PROTECT(rn = allocVector(INTSXP, 2)); /* 33 */
    INTEGER(rn)[0] = NA_INTEGER; INTEGER(rn)[1] = -n_tdcov_rows;
    setAttrib(tdc_df, R_RowNamesSymbol, rn);
    UNPROTECT(1); /* rn; stack = 32 */

    /* Build result list — 4 elements (33, 34) */
    PROTECT(result = allocVector(VECSXP, 4)); /* 33 */
    PROTECT(res_nm = allocVector(STRSXP, 4)); /* 34 */
    SET_VECTOR_ELT(result, 0, pop_df);  SET_STRING_ELT(res_nm, 0, mkChar("pop"));
    SET_VECTOR_ELT(result, 1, tic_df);  SET_STRING_ELT(res_nm, 1, mkChar("time_ind_covariate"));
    SET_VECTOR_ELT(result, 2, com_df);  SET_STRING_ELT(res_nm, 2, mkChar("community"));
    SET_VECTOR_ELT(result, 3, tdc_df);  SET_STRING_ELT(res_nm, 3, mkChar("time_dep_covariate"));
    setAttrib(result, R_NamesSymbol, res_nm);

    UNPROTECT(34);
    return result;
}

/* ---------------------------------------------------------
 * r_transtat
 *
 * Run the full TranStat estimation pipeline.
 * Reads config.file from the given path, then runs maximum
 * likelihood estimation; results are written to the output
 * directory specified in config.file.
 *
 * Arguments:
 *   r_config_file   – full path to config.file
 *   r_serial_number – integer identifying this run (used in output file names)
 *   r_seed          – integer RNG seed for Monte Carlo procedures
 *
 * Returns NULL invisibly; results are written to files.
 * --------------------------------------------------------- */
SEXP r_transtat(SEXP r_config_file, SEXP r_serial_number, SEXP r_seed)
{
    int i, j, k, l;
    int len_inc, len_inf;
    FILE *file;
    char file_name[512];

    serial_number = INTEGER(r_serial_number)[0];
    seed          = (long)INTEGER(r_seed)[0];
    mt_init(seed);

    if ((file = fopen(CHAR(STRING_ELT(r_config_file, 0)), "r")) == NULL)
        error("Cannot open config file: %s", CHAR(STRING_ELT(r_config_file, 0)));
    cfg_pars = get_cfg_pars(file);
    fclose(file);

    for (i = 0; i < cfg_pars.n_inc; i++) {
        cfg_pars.max_incubation  = cfg_pars.max_inc[i];
        cfg_pars.min_incubation  = cfg_pars.min_inc[i];
        len_inc = cfg_pars.max_inc[i] - cfg_pars.min_inc[i] + 1;
        cfg_pars.prob_incubation = (double *)malloc((size_t)(len_inc * sizeof(double)));
        for (k = 0; k < len_inc; k++)
            cfg_pars.prob_incubation[k] = cfg_pars.prob_inc[i][k];

        for (j = 0; j < cfg_pars.n_inf; j++) {
            cfg_pars.upper_infectious = cfg_pars.upper_inf[j];
            cfg_pars.lower_infectious = cfg_pars.lower_inf[j];
            len_inf = cfg_pars.upper_inf[j] - cfg_pars.lower_inf[j] + 1;
            cfg_pars.prob_infectious  = (double *)malloc((size_t)(len_inf * sizeof(double)));
            for (l = 0; l < len_inf; l++)
                cfg_pars.prob_infectious[l] = cfg_pars.prob_inf[j][l];

            if (cfg_pars.simplify_output == 0) {
                sprintf(file_name, "%soutput.txt", cfg_pars.path_out);
                if ((file = fopen(file_name, "a")) == NULL)
                    file = fopen(file_name, "w");
                if (file) {
                    fprintf(file, "\n\n Serial number=%d\n", serial_number);
                    fprintf(file,
                            "=====================================================================\n");
                    fprintf(file, "%d: min_inc=%d  max_inc=%d\n",
                            i+1, cfg_pars.min_incubation, cfg_pars.max_incubation);
                    fprintf(file, "%d: lower_inf=%d  upper_inf=%d\n",
                            j+1, cfg_pars.lower_infectious, cfg_pars.upper_infectious);
                    fprintf(file,
                            "=====================================================================\n");
                    fclose(file);
                }
            }

            core(i, j, 0);
            free(cfg_pars.prob_infectious);
        }
        free(cfg_pars.prob_incubation);
    }

    /* Free cfg_pars memory (mirrors main.c cleanup) */
    free(cfg_pars.min_inc);
    free(cfg_pars.max_inc);
    free(cfg_pars.lower_inf);
    free(cfg_pars.upper_inf);
    for (i = 0; i < cfg_pars.n_inc; i++) free(cfg_pars.prob_inc[i]);
    free(cfg_pars.prob_inc);
    for (i = 0; i < cfg_pars.n_inf; i++) free(cfg_pars.prob_inf[i]);
    free(cfg_pars.prob_inf);
    free(cfg_pars.c2p_covariate);
    free(cfg_pars.sus_p2p_covariate);
    free(cfg_pars.inf_p2p_covariate);
    free(cfg_pars.pat_covariate);
    free(cfg_pars.imm_covariate);
    free_2d_array_int(cfg_pars.interaction);
    free_2d_array_double(cfg_pars.ini_par_effective);
    free(cfg_pars.lower_search_bound);
    free(cfg_pars.upper_search_bound);
    free(cfg_pars.converge_criteria);
    free_2d_array_double(cfg_pars.SAR_sus_time_ind_covariate);
    free_2d_array_double(cfg_pars.SAR_inf_time_ind_covariate);
    free_3d_array_double(cfg_pars.SAR_sus_time_dep_covariate);
    free_3d_array_double(cfg_pars.SAR_inf_time_dep_covariate);
    free(cfg_pars.par_fixed_id);
    free(cfg_pars.par_fixed_value);
    free(cfg_pars.sim_par_effective);
    free(cfg_pars.effective_lower_infectious);
    free(cfg_pars.effective_upper_infectious);
    if (cfg_pars.par_equiclass != NULL) {
        for (i = 0; i < cfg_pars.n_par_equiclass; i++)
            free(cfg_pars.par_equiclass[i].member);
        free(cfg_pars.par_equiclass);
    }
    free(cfg_pars.c2p_group);

    return R_NilValue;
}

/* =========================================================
 * Helpers for extracting named elements from R lists
 * ========================================================= */

static SEXP cb_get(SEXP lst, const char *nm) {
    SEXP names = getAttrib(lst, R_NamesSymbol);
    int n = length(names);
    for (int k = 0; k < n; k++)
        if (strcmp(CHAR(STRING_ELT(names, k)), nm) == 0)
            return VECTOR_ELT(lst, k);
    return R_NilValue;
}
static int cb_gi(SEXP lst, const char *nm) {
    SEXP e = cb_get(lst, nm);
    if (e == R_NilValue || length(e) == 0) return 0;
    return INTEGER(e)[0];
}
static double cb_gd(SEXP lst, const char *nm) {
    SEXP e = cb_get(lst, nm);
    if (e == R_NilValue || length(e) == 0) return 0.0;
    if (isReal(e))    return REAL(e)[0];
    if (isInteger(e)) return (double)INTEGER(e)[0];
    return 0.0;
}

/* ---------------------------------------------------------
 * r_simulate_single
 *
 * Run one epidemic simulation and return the updated pop_list
 * plus summary statistics.  Unlike r_transtat(), this function
 * populates all C global state directly from R data frames
 * (no temporary files), so the post-simulation time-dependent
 * covariate matrix (e.g. treatment status) is available to R.
 *
 * Arguments:
 *   r_pop       – pop data frame (13 cols, see gen_population)
 *   r_tic       – time_ind_covariate data frame (id, value)
 *   r_community – community data frame (5 cols)
 *   r_tdc       – time_dep_covariate data frame (id, day_start, day_stop, value)
 *   r_cfg       – simulation config list (output of read_config, with
 *                 sim_par_effective already logit/log-transformed)
 *   r_i_inc     – 1-based index into incubation-period group
 *   r_i_inf     – 1-based index into infectious-period group
 *   r_seed      – integer RNG seed
 *
 * Returns a named list:
 *   $pop, $time_ind_covariate, $community, $time_dep_covariate
 *   (updated data frames)
 *   $n_index, $n_secondary_inf, $n_secondary_sym, $n_secondary_asym,
 *   $n_escaped, $n_preimmune  (summary counts)
 * --------------------------------------------------------- */
SEXP r_simulate_single(SEXP r_pop, SEXP r_tic, SEXP r_community, SEXP r_tdc,
                       SEXP r_cfg, SEXP r_i_inc, SEXP r_i_inf, SEXP r_seed)
{
    int h, i, j, k, m, r, t;
    int i_inc, i_inf;
    int loc_p_size, loc_n_community, n_tdc_rows;
    int len_inc, len_inf;
    int n_sym_idx, n_asym_idx, n_sym_sec, n_asym_sec;
    int n_esc_attacked, n_imm;
    CONTACT *ptr_contact, *ptr2_contact;
    int n_protect = 0;

    /* --- 1. Extract indices and seed --- */
    i_inc = INTEGER(r_i_inc)[0] - 1;  /* 0-based */
    i_inf = INTEGER(r_i_inf)[0] - 1;
    seed  = (long)INTEGER(r_seed)[0];
    mt_init(seed);

    loc_p_size       = length(VECTOR_ELT(r_pop, 0));
    loc_n_community  = length(VECTOR_ELT(r_community, 0));
    n_tdc_rows       = length(VECTOR_ELT(r_tdc, 0));

    /* --- 2. Build cfg_pars from r_cfg --- */
    memset(&cfg_pars, 0, sizeof(CFG_PARS));

    cfg_pars.n_b_mode             = cb_gi(r_cfg, "n_b_mode");
    cfg_pars.n_p_mode             = cb_gi(r_cfg, "n_p_mode");
    cfg_pars.n_u_mode             = cb_gi(r_cfg, "n_u_mode");
    cfg_pars.n_q_mode             = cb_gi(r_cfg, "n_q_mode");
    cfg_pars.n_time_ind_covariate = cb_gi(r_cfg, "n_time_ind_covariate");
    cfg_pars.n_time_dep_covariate = cb_gi(r_cfg, "n_time_dep_covariate");
    cfg_pars.n_covariate          = cfg_pars.n_time_ind_covariate + cfg_pars.n_time_dep_covariate;
    cfg_pars.n_c2p_covariate      = cb_gi(r_cfg, "n_c2p_covariate");
    cfg_pars.n_sus_p2p_covariate  = cb_gi(r_cfg, "n_sus_p2p_covariate");
    cfg_pars.n_inf_p2p_covariate  = cb_gi(r_cfg, "n_inf_p2p_covariate");
    cfg_pars.n_int_p2p_covariate  = cb_gi(r_cfg, "n_int_p2p_covariate");
    cfg_pars.n_p2p_covariate      = cb_gi(r_cfg, "n_p2p_covariate");
    cfg_pars.n_pat_covariate      = cb_gi(r_cfg, "n_pat_covariate");
    cfg_pars.n_imm_covariate      = cb_gi(r_cfg, "n_imm_covariate");
    cfg_pars.n_par                = cb_gi(r_cfg, "n_par");
    cfg_pars.n_par_equiclass      = cb_gi(r_cfg, "n_par_equiclass");
    cfg_pars.n_par_fixed          = cb_gi(r_cfg, "n_par_fixed");

    /* incubation period: select i_inc-th group */
    {
        SEXP min_inc_r  = cb_get(r_cfg, "min_inc");
        SEXP max_inc_r  = cb_get(r_cfg, "max_inc");
        SEXP prob_inc_r = cb_get(r_cfg, "prob_inc");   /* R list */
        cfg_pars.min_incubation = INTEGER(min_inc_r)[i_inc];
        cfg_pars.max_incubation = INTEGER(max_inc_r)[i_inc];
        len_inc = cfg_pars.max_incubation - cfg_pars.min_incubation + 1;
        cfg_pars.prob_incubation = (double *)malloc((size_t)(len_inc * sizeof(double)));
        SEXP pvec = VECTOR_ELT(prob_inc_r, i_inc);
        for (k = 0; k < len_inc; k++) cfg_pars.prob_incubation[k] = REAL(pvec)[k];
    }

    /* infectious period: select i_inf-th group, then modify for pre-onset */
    {
        SEXP lower_inf_r  = cb_get(r_cfg, "lower_inf");
        SEXP upper_inf_r  = cb_get(r_cfg, "upper_inf");
        SEXP prob_inf_r   = cb_get(r_cfg, "prob_inf");   /* R list */
        cfg_pars.lower_infectious = INTEGER(lower_inf_r)[i_inf];
        cfg_pars.upper_infectious = INTEGER(upper_inf_r)[i_inf];
        len_inf = cfg_pars.upper_infectious - cfg_pars.lower_infectious + 1;
        cfg_pars.prob_infectious = (double *)malloc((size_t)(len_inf * sizeof(double)));
        SEXP pvec = VECTOR_ELT(prob_inf_r, i_inf);
        for (k = 0; k < len_inf; k++) cfg_pars.prob_infectious[k] = REAL(pvec)[k];

        /* Adjust prob_infectious for pre-symptom-onset days (lower_infectious < 0) */
        if (cfg_pars.lower_infectious < 0) {
            int l2 = cfg_pars.max_incubation - cfg_pars.min_incubation + 1;
            double *pdf_inc = (double *)malloc((size_t)(l2 * sizeof(double)));
            double *sdf_inc = (double *)malloc((size_t)(l2 * sizeof(double)));
            for (k = 0; k < l2; k++) pdf_inc[k] = cfg_pars.prob_incubation[k];
            sdf_inc[0] = 1.0 - pdf_inc[0];
            for (k = 1; k < l2; k++)
                sdf_inc[k] = (k == l2-1) ? 0.0 : sdf_inc[k-1] - pdf_inc[k];
            for (t = cfg_pars.lower_infectious; t < 0; t++) {
                r = t - cfg_pars.lower_infectious;
                int linc = -t;
                if (linc >= cfg_pars.min_incubation) {
                    if (linc <= cfg_pars.max_incubation)
                        cfg_pars.prob_infectious[r] *= sdf_inc[linc - cfg_pars.min_incubation]
                                                     + pdf_inc[linc - cfg_pars.min_incubation];
                    else
                        cfg_pars.prob_infectious[r] = 0.0;
                }
            }
            free(pdf_inc); free(sdf_inc);
        }
    }

    /* covariate index arrays */
    {
        SEXP e;
        int n;
        n = cfg_pars.n_c2p_covariate;
        if (n > 0) {
            cfg_pars.c2p_covariate = (int *)malloc((size_t)(n * sizeof(int)));
            e = cb_get(r_cfg, "c2p_covariate");
            for (k = 0; k < n; k++) cfg_pars.c2p_covariate[k] = INTEGER(e)[k];
        }
        n = cfg_pars.n_sus_p2p_covariate;
        if (n > 0) {
            cfg_pars.sus_p2p_covariate = (int *)malloc((size_t)(n * sizeof(int)));
            e = cb_get(r_cfg, "sus_p2p_covariate");
            for (k = 0; k < n; k++) cfg_pars.sus_p2p_covariate[k] = INTEGER(e)[k];
        }
        n = cfg_pars.n_inf_p2p_covariate;
        if (n > 0) {
            cfg_pars.inf_p2p_covariate = (int *)malloc((size_t)(n * sizeof(int)));
            e = cb_get(r_cfg, "inf_p2p_covariate");
            for (k = 0; k < n; k++) cfg_pars.inf_p2p_covariate[k] = INTEGER(e)[k];
        }
        n = cfg_pars.n_int_p2p_covariate;
        if (n > 0) {
            make_2d_array_int(&cfg_pars.interaction, n, 2, 0);
            SEXP rint = cb_get(r_cfg, "interaction");   /* list of c(sus,inf) pairs */
            for (k = 0; k < n; k++) {
                SEXP pair = VECTOR_ELT(rint, k);
                cfg_pars.interaction[k][0] = INTEGER(pair)[0];
                cfg_pars.interaction[k][1] = INTEGER(pair)[1];
            }
        }
        n = cfg_pars.n_pat_covariate;
        if (n > 0) {
            cfg_pars.pat_covariate = (int *)malloc((size_t)(n * sizeof(int)));
            e = cb_get(r_cfg, "pat_covariate");
            for (k = 0; k < n; k++) cfg_pars.pat_covariate[k] = INTEGER(e)[k];
        }
        n = cfg_pars.n_imm_covariate;
        if (n > 0) {
            cfg_pars.imm_covariate = (int *)malloc((size_t)(n * sizeof(int)));
            e = cb_get(r_cfg, "imm_covariate");
            for (k = 0; k < n; k++) cfg_pars.imm_covariate[k] = INTEGER(e)[k];
        }
    }

    /* parameter equiclasses */
    {
        int n_peq = cfg_pars.n_par_equiclass;
        cfg_pars.par_equiclass = NULL;
        if (n_peq > 0) {
            cfg_pars.par_equiclass = (PAR_EQUICLASS *)malloc((size_t)(n_peq * sizeof(PAR_EQUICLASS)));
            SEXP r_peq = cb_get(r_cfg, "par_equiclass");
            for (k = 0; k < n_peq; k++) {
                SEXP cls  = VECTOR_ELT(r_peq, k);
                int  sz   = cb_gi(cls, "size");
                cfg_pars.par_equiclass[k].size   = sz;
                cfg_pars.par_equiclass[k].member = (int *)malloc((size_t)(sz * sizeof(int)));
                SEXP memb = cb_get(cls, "member");
                for (j = 0; j < sz; j++)
                    cfg_pars.par_equiclass[k].member[j] = INTEGER(memb)[j];
            }
        }
    }

    /* sim_par_effective (already logit/log-transformed by R wrapper) */
    {
        int n_peq = cfg_pars.n_par_equiclass;
        cfg_pars.sim_par_effective = (double *)malloc((size_t)(n_peq * sizeof(double)));
        SEXP spe = cb_get(r_cfg, "sim_par_effective_tr");
        for (k = 0; k < n_peq; k++) cfg_pars.sim_par_effective[k] = REAL(spe)[k];
    }

    /* fixed parameters */
    cfg_pars.par_fixed_id    = NULL;
    cfg_pars.par_fixed_value = NULL;
    if (cfg_pars.n_par_fixed > 0) {
        cfg_pars.par_fixed_id    = (int *)malloc((size_t)(cfg_pars.n_par_fixed * sizeof(int)));
        cfg_pars.par_fixed_value = (double *)malloc((size_t)(cfg_pars.n_par_fixed * sizeof(double)));
        SEXP fid = cb_get(r_cfg, "par_fixed_id");
        SEXP fvl = cb_get(r_cfg, "par_fixed_value_tr");
        for (k = 0; k < cfg_pars.n_par_fixed; k++) {
            cfg_pars.par_fixed_id[k]    = INTEGER(fid)[k];
            cfg_pars.par_fixed_value[k] = REAL(fvl)[k];
        }
    }

    /* symptom-induced covariate flags */
    cfg_pars.PreIllness_as_covariate = cb_gi(r_cfg, "PreIllness_as_covariate");
    cfg_pars.PreIllness_covariate_id = cb_gi(r_cfg, "PreIllness_covariate_id");
    cfg_pars.Illness_as_covariate    = cb_gi(r_cfg, "Illness_as_covariate");
    cfg_pars.Illness_covariate_id    = cb_gi(r_cfg, "Illness_covariate_id");
    cfg_pars.RxIllness_as_covariate  = cb_gi(r_cfg, "RxIllness_as_covariate");
    cfg_pars.RxIllness_covariate_id  = cb_gi(r_cfg, "RxIllness_covariate_id");
    cfg_pars.RxIllness_prob          = cb_gd(r_cfg, "RxIllness_prob");
    cfg_pars.RxIllness_duration      = cb_gi(r_cfg, "RxIllness_duration");
    cfg_pars.RxIllness_index_only    = cb_gi(r_cfg, "RxIllness_index_only");

    /* contact generation flags */
    cfg_pars.generate_c2p_contact                = cb_gi(r_cfg, "generate_c2p_contact");
    cfg_pars.generate_p2p_contact                = cb_gi(r_cfg, "generate_p2p_contact");
    cfg_pars.common_contact_history_within_community = cb_gi(r_cfg, "common_contact_history_within_community");
    cfg_pars.c2p_offset                          = cb_gi(r_cfg, "c2p_offset");
    cfg_pars.p2p_offset                          = cb_gi(r_cfg, "p2p_offset");

    /* simulation control */
    cfg_pars.simulation          = 1;
    cfg_pars.n_simulation        = 1;
    cfg_pars.simulation_only     = 1;
    cfg_pars.EM                  = 0;
    cfg_pars.preset_index        = cb_gi(r_cfg, "preset_index");
    cfg_pars.adjust_for_left_truncation  = cb_gi(r_cfg, "adjust_for_left_truncation");
    cfg_pars.adjust_for_right_censoring = cb_gi(r_cfg, "adjust_for_right_censoring");
    cfg_pars.idx_initiated_followup_in_simulation      = cb_gi(r_cfg, "idx_initiated_followup_in_simulation");
    cfg_pars.followup_duration_after_idx_in_simulation = cb_gi(r_cfg, "followup_duration_after_idx_in_simulation");
    cfg_pars.asym_effect_sim     = cb_gd(r_cfg, "asym_effect_sim");
    cfg_pars.silent_run          = cb_gi(r_cfg, "silent_run");

    /* unused-for-simulation fields: NULL/0 */
    cfg_pars.SAR_n_covariate_sets       = 0;
    cfg_pars.SAR_sus_time_ind_covariate = NULL;
    cfg_pars.SAR_inf_time_ind_covariate = NULL;
    cfg_pars.SAR_sus_time_dep_covariate = NULL;
    cfg_pars.SAR_inf_time_dep_covariate = NULL;
    cfg_pars.effective_lower_infectious = NULL;
    cfg_pars.effective_upper_infectious = NULL;
    cfg_pars.c2p_group                  = NULL;
    cfg_pars.n_c2p_group                = 0;
    cfg_pars.use_index_cases_to_improve_b = 0;
    cfg_pars.n_inc                      = 1;
    cfg_pars.n_inf                      = 1;

    /* --- 3. Set global sizes --- */
    p_size      = loc_p_size;
    n_community = loc_n_community;
    importance_weight  = NULL;
    pdf_incubation = sdf_incubation = NULL;
    size_sample_states = 0;
    max_epi_duration   = 0;

    /* --- 4. Allocate and fill people[] --- */
    people = (PEOPLE *)malloc((size_t)(p_size * sizeof(PEOPLE)));
    if (!people) error("Cannot allocate people array in r_simulate_single");

    {
        /* Column pointers into r_pop */
        int    *p_id   = INTEGER(VECTOR_ELT(r_pop,  0));
        int    *p_comm = INTEGER(VECTOR_ELT(r_pop,  1));
        int    *p_pimm = INTEGER(VECTOR_ELT(r_pop,  2));
        int    *p_inf  = INTEGER(VECTOR_ELT(r_pop,  3));
        int    *p_sym  = INTEGER(VECTOR_ELT(r_pop,  4));
        int    *p_dill = INTEGER(VECTOR_ELT(r_pop,  5));
        int    *p_exit = INTEGER(VECTOR_ELT(r_pop,  6));
        int    *p_dex  = INTEGER(VECTOR_ELT(r_pop,  7));
        int    *p_idx  = INTEGER(VECTOR_ELT(r_pop,  8));
        int    *p_umod = INTEGER(VECTOR_ELT(r_pop,  9));
        int    *p_qmod = INTEGER(VECTOR_ELT(r_pop, 10));
        double *p_wt   = REAL   (VECTOR_ELT(r_pop, 11));
        int    *p_ign  = INTEGER(VECTOR_ELT(r_pop, 12));

        for (i = 0; i < p_size; i++) {
            people[i].id         = p_id[i];
            people[i].community  = p_comm[i];
            people[i].pre_immune = p_pimm[i];
            people[i].infection  = p_inf[i];
            people[i].symptom    = p_sym[i];
            people[i].day_ill    = p_dill[i];
            people[i].exit       = p_exit[i];
            people[i].day_exit   = p_dex[i];
            people[i].idx        = p_idx[i];
            people[i].u_mode     = p_umod[i];
            people[i].q_mode     = p_qmod[i];
            people[i].weight     = p_wt[i];
            people[i].ignore     = p_ign[i];
            /* derived fields */
            people[i].day_infection       = CB_MISSING;
            people[i].day_infection_lower = CB_MISSING;
            people[i].day_infection_upper = CB_MISSING;
            people[i].day_infective_lower = CB_MISSING;
            people[i].day_infective_upper = CB_MISSING;
            people[i].final_risk_day      = CB_MISSING;
            people[i].any_exposure        = 0;
            people[i].current_state       = 0;
            people[i].size_possible_states = 0;
            /* pointers: allocated below */
            people[i].time_ind_covariate = NULL;
            people[i].time_dep_covariate = NULL;
            people[i].contact_history    = NULL;
            people[i].risk_history       = NULL;
            people[i].risk_class         = NULL;
            people[i].imm_covariate      = NULL;
            people[i].pat_covariate      = NULL;
            people[i].possible_states    = NULL;
        }
    }

    /* --- 5. Allocate community[] --- */
    community = (COMMUNITY *)malloc((size_t)(n_community * sizeof(COMMUNITY)));
    if (!community) { free(people); error("Cannot allocate community array"); }

    {
        int *c_id  = INTEGER(VECTOR_ELT(r_community, 0));
        int *c_ds  = INTEGER(VECTOR_ELT(r_community, 1));
        int *c_dp  = INTEGER(VECTOR_ELT(r_community, 2));
        int *c_dfu = INTEGER(VECTOR_ELT(r_community, 3));
        int *c_c2p = INTEGER(VECTOR_ELT(r_community, 4));

        for (h = 0; h < n_community; h++) {
            community[h].id                  = c_id[h];
            community[h].day_epi_start       = c_ds[h];
            community[h].day_epi_stop        = c_dp[h];
            community[h].day_last_followup   = c_dfu[h];
            community[h].c2p_group           = c_c2p[h];
            community[h].epi_duration        = c_dp[h] - c_ds[h] + 1;
            if (community[h].epi_duration > max_epi_duration)
                max_epi_duration = community[h].epi_duration;
            community[h].size                = 0;
            community[h].size_idx            = 0;
            community[h].earliest_idx_day_ill = CB_MISSING;
            community[h].latest_idx_day_ill   = CB_MISSING;
            community[h].counter             = 0;
            community[h].ignore              = 0;
            community[h].member              = NULL;
            community[h].idx                 = NULL;
            community[h].member_impute       = NULL;
            community[h].size_impute         = 0;
            community[h].size_possible_states = 1;
            community[h].contact_history     = NULL;
            community[h].risk_class          = NULL;
            community[h].risk_class_rear     = NULL;
            community[h].sample_states       = NULL;
            community[h].sample_states_rear  = NULL;
            community[h].list_states         = NULL;
            community[h].list_states_rear    = NULL;
        }
    }

    /* --- 6. Build community membership --- */
    for (i = 0; i < p_size; i++) {
        if (people[i].ignore == 0) {
            h = people[i].community;
            community[h].size++;
        }
    }
    for (h = 0; h < n_community; h++) {
        if (community[h].size > 0)
            community[h].member = (int *)malloc((size_t)(community[h].size * sizeof(int)));
    }
    for (h = 0; h < n_community; h++) community[h].counter = 0;
    for (i = 0; i < p_size; i++) {
        if (people[i].ignore == 0) {
            h = people[i].community;
            community[h].member[community[h].counter++] = i;
        }
    }

    /* --- 7. Allocate per-person covariate arrays and fill --- */
    {
        int n_tic = cfg_pars.n_time_ind_covariate;
        int n_tdc = cfg_pars.n_time_dep_covariate;
        for (i = 0; i < p_size; i++) {
            if (people[i].ignore == 0) {
                h = people[i].community;
                if (n_tic > 0) make_1d_array_double(&people[i].time_ind_covariate, n_tic, 0.0);
                if (n_tdc > 0) make_2d_array_double(&people[i].time_dep_covariate,
                                                    community[h].epi_duration, n_tdc, 0.0);
            }
        }
    }

    /* fill time-independent covariates from r_tic */
    {
        int   *tic_id  = INTEGER(VECTOR_ELT(r_tic, 0));
        double *tic_v  = REAL   (VECTOR_ELT(r_tic, 1));
        int n_tic_rows = length(VECTOR_ELT(r_tic, 0));
        /* r_tic has one row per person with all n_time_ind_covariate values.
           Currently gen_population() stores a single covariate per row. */
        for (m = 0; m < n_tic_rows; m++) {
            i = tic_id[m];
            if (i >= 0 && i < p_size && people[i].ignore == 0 &&
                people[i].time_ind_covariate != NULL)
                people[i].time_ind_covariate[0] = tic_v[m];
        }
    }

    /* fill time-dependent covariates from r_tdc */
    {
        int    *tdc_id  = INTEGER(VECTOR_ELT(r_tdc, 0));
        int    *tdc_ds  = INTEGER(VECTOR_ELT(r_tdc, 1));
        int    *tdc_dp  = INTEGER(VECTOR_ELT(r_tdc, 2));
        double *tdc_v   = REAL   (VECTOR_ELT(r_tdc, 3));
        for (m = 0; m < n_tdc_rows; m++) {
            i = tdc_id[m];
            if (i < 0 || i >= p_size || people[i].ignore != 0 ||
                people[i].time_dep_covariate == NULL) continue;
            h = people[i].community;
            for (t = tdc_ds[m]; t <= tdc_dp[m]; t++) {
                if (t >= community[h].day_epi_start && t <= community[h].day_epi_stop) {
                    r = t - community[h].day_epi_start;
                    people[i].time_dep_covariate[r][0] = tdc_v[m];
                }
            }
        }
    }

    /* set day_exit for people with exit==0 (no censoring → very large day) */
    for (i = 0; i < p_size; i++) {
        if (people[i].exit == 0) people[i].day_exit = 1000000; /* INFINITY_INTEGER */
        if (people[i].infection == 0) people[i].day_ill = CB_MISSING;
    }

    /* --- 8. Allocate community contact histories (common_contact_history_within_community=1) --- */
    for (h = 0; h < n_community; h++) {
        if (community[h].size > 0 && cfg_pars.common_contact_history_within_community == 1) {
            community[h].contact_history = (CONTACT_HISTORY *)malloc(
                (size_t)(community[h].epi_duration * sizeof(CONTACT_HISTORY)));
            for (t = community[h].day_epi_start; t <= community[h].day_epi_stop; t++) {
                r = t - community[h].day_epi_start;
                community[h].contact_history[r].c2p_contact      = NULL;
                community[h].contact_history[r].c2p_contact_rear = NULL;
                community[h].contact_history[r].p2p_contact      = NULL;
                community[h].contact_history[r].p2p_contact_rear = NULL;
            }
        }
    }

    /* --- 9. Allocate global parameter arrays --- */
    create_arrays();

    /* --- 10. Auto-generate contact histories --- */
    if (cfg_pars.generate_c2p_contact == 1) {
        for (h = 0; h < n_community; h++) {
            if (community[h].size > 0)
                add_c2p_contact_history_to_community(h, community[h].day_epi_start,
                                                     community[h].day_epi_stop, 0, 0.0);
        }
    }
    if (cfg_pars.generate_p2p_contact == 1) {
        for (h = 0; h < n_community; h++) {
            if (community[h].size > 0)
                add_p2p_contact_history_to_community(h, community[h].day_epi_start,
                                                     community[h].day_epi_stop, 0, 0.0);
        }
    }

    /* --- 11. Run one simulation --- */
    simulate();

    /* --- 12. Summary statistics --- */
    n_sym_idx = n_asym_idx = n_sym_sec = n_asym_sec = n_esc_attacked = n_imm = 0;
    for (i = 0; i < p_size; i++) {
        h = people[i].community;
        if (people[i].pre_immune == 1) { n_imm++; continue; }
        if (people[i].infection == 0) {
            if (community[h].size_idx > 0) n_esc_attacked++;
            continue;
        }
        if (people[i].idx == 1) {
            if (people[i].symptom == 1) n_sym_idx++; else n_asym_idx++;
        } else {
            if (people[i].symptom == 1) n_sym_sec++; else n_asym_sec++;
        }
    }

    /* --- 13. Build output data frames --- */

    /* pop: 14 columns (original 13 + day_infection) */
    SEXP o_pop_id, o_pop_comm, o_pop_pimm, o_pop_inf, o_pop_sym;
    SEXP o_pop_dill, o_pop_exit, o_pop_dex, o_pop_idx, o_pop_umod, o_pop_qmod;
    SEXP o_pop_wt, o_pop_ign, o_pop_dinf;
    PROTECT(o_pop_id   = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_pop_comm = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_pop_pimm = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_pop_inf  = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_pop_sym  = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_pop_dill = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_pop_exit = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_pop_dex  = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_pop_idx  = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_pop_umod = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_pop_qmod = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_pop_wt   = allocVector(REALSXP, p_size)); n_protect++;
    PROTECT(o_pop_ign  = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_pop_dinf = allocVector(INTSXP,  p_size)); n_protect++;

    for (i = 0; i < p_size; i++) {
        INTEGER(o_pop_id)[i]   = people[i].id;
        INTEGER(o_pop_comm)[i] = people[i].community;
        INTEGER(o_pop_pimm)[i] = people[i].pre_immune;
        INTEGER(o_pop_inf)[i]  = people[i].infection;
        INTEGER(o_pop_sym)[i]  = people[i].symptom;
        INTEGER(o_pop_dill)[i] = (people[i].day_ill == CB_MISSING) ? NA_INTEGER : people[i].day_ill;
        INTEGER(o_pop_exit)[i] = people[i].exit;
        INTEGER(o_pop_dex)[i]  = (people[i].day_exit == 1000000) ? CB_MISSING : people[i].day_exit;
        INTEGER(o_pop_idx)[i]  = people[i].idx;
        INTEGER(o_pop_umod)[i] = people[i].u_mode;
        INTEGER(o_pop_qmod)[i] = people[i].q_mode;
        REAL(o_pop_wt)[i]      = people[i].weight;
        INTEGER(o_pop_ign)[i]  = people[i].ignore;
        INTEGER(o_pop_dinf)[i] = (people[i].day_infection == CB_MISSING) ? NA_INTEGER : people[i].day_infection;
    }

    SEXP pop_df, pop_nm, pop_rn;
    PROTECT(pop_df = allocVector(VECSXP, 14)); n_protect++;
    PROTECT(pop_nm = allocVector(STRSXP, 14)); n_protect++;
    SET_VECTOR_ELT(pop_df,  0, o_pop_id);   SET_STRING_ELT(pop_nm,  0, mkChar("id"));
    SET_VECTOR_ELT(pop_df,  1, o_pop_comm); SET_STRING_ELT(pop_nm,  1, mkChar("community"));
    SET_VECTOR_ELT(pop_df,  2, o_pop_pimm); SET_STRING_ELT(pop_nm,  2, mkChar("pre_immune"));
    SET_VECTOR_ELT(pop_df,  3, o_pop_inf);  SET_STRING_ELT(pop_nm,  3, mkChar("infection"));
    SET_VECTOR_ELT(pop_df,  4, o_pop_sym);  SET_STRING_ELT(pop_nm,  4, mkChar("symptom"));
    SET_VECTOR_ELT(pop_df,  5, o_pop_dill); SET_STRING_ELT(pop_nm,  5, mkChar("day_ill"));
    SET_VECTOR_ELT(pop_df,  6, o_pop_exit); SET_STRING_ELT(pop_nm,  6, mkChar("exit"));
    SET_VECTOR_ELT(pop_df,  7, o_pop_dex);  SET_STRING_ELT(pop_nm,  7, mkChar("day_exit"));
    SET_VECTOR_ELT(pop_df,  8, o_pop_idx);  SET_STRING_ELT(pop_nm,  8, mkChar("idx"));
    SET_VECTOR_ELT(pop_df,  9, o_pop_umod); SET_STRING_ELT(pop_nm,  9, mkChar("u_mode"));
    SET_VECTOR_ELT(pop_df, 10, o_pop_qmod); SET_STRING_ELT(pop_nm, 10, mkChar("q_mode"));
    SET_VECTOR_ELT(pop_df, 11, o_pop_wt);   SET_STRING_ELT(pop_nm, 11, mkChar("weight"));
    SET_VECTOR_ELT(pop_df, 12, o_pop_ign);  SET_STRING_ELT(pop_nm, 12, mkChar("ignore"));
    SET_VECTOR_ELT(pop_df, 13, o_pop_dinf); SET_STRING_ELT(pop_nm, 13, mkChar("day_infection"));
    setAttrib(pop_df, R_NamesSymbol, pop_nm);
    setAttrib(pop_df, R_ClassSymbol, mkString("data.frame"));
    PROTECT(pop_rn = allocVector(INTSXP, 2)); n_protect++;
    INTEGER(pop_rn)[0] = NA_INTEGER; INTEGER(pop_rn)[1] = -p_size;
    setAttrib(pop_df, R_RowNamesSymbol, pop_rn);

    /* time_ind_covariate: 2 columns (id, value) */
    SEXP o_tic_id, o_tic_v, tic_df, tic_nm, tic_rn;
    PROTECT(o_tic_id = allocVector(INTSXP,  p_size)); n_protect++;
    PROTECT(o_tic_v  = allocVector(REALSXP, p_size)); n_protect++;
    for (i = 0; i < p_size; i++) {
        INTEGER(o_tic_id)[i] = people[i].id;
        REAL(o_tic_v)[i] = (people[i].time_ind_covariate != NULL) ?
                            people[i].time_ind_covariate[0] : 0.0;
    }
    PROTECT(tic_df = allocVector(VECSXP, 2)); n_protect++;
    PROTECT(tic_nm = allocVector(STRSXP, 2)); n_protect++;
    SET_VECTOR_ELT(tic_df, 0, o_tic_id); SET_STRING_ELT(tic_nm, 0, mkChar("id"));
    SET_VECTOR_ELT(tic_df, 1, o_tic_v);  SET_STRING_ELT(tic_nm, 1, mkChar("value"));
    setAttrib(tic_df, R_NamesSymbol, tic_nm);
    setAttrib(tic_df, R_ClassSymbol, mkString("data.frame"));
    PROTECT(tic_rn = allocVector(INTSXP, 2)); n_protect++;
    INTEGER(tic_rn)[0] = NA_INTEGER; INTEGER(tic_rn)[1] = -p_size;
    setAttrib(tic_df, R_RowNamesSymbol, tic_rn);

    /* community: 7 columns */
    SEXP o_com_id, o_com_ds, o_com_dp, o_com_dfu, o_com_c2p, o_com_eidx, o_com_ign;
    PROTECT(o_com_id   = allocVector(INTSXP, n_community)); n_protect++;
    PROTECT(o_com_ds   = allocVector(INTSXP, n_community)); n_protect++;
    PROTECT(o_com_dp   = allocVector(INTSXP, n_community)); n_protect++;
    PROTECT(o_com_dfu  = allocVector(INTSXP, n_community)); n_protect++;
    PROTECT(o_com_c2p  = allocVector(INTSXP, n_community)); n_protect++;
    PROTECT(o_com_eidx = allocVector(INTSXP, n_community)); n_protect++;
    PROTECT(o_com_ign  = allocVector(INTSXP, n_community)); n_protect++;
    for (h = 0; h < n_community; h++) {
        INTEGER(o_com_id)[h]   = community[h].id;
        INTEGER(o_com_ds)[h]   = community[h].day_epi_start;
        INTEGER(o_com_dp)[h]   = community[h].day_epi_stop;
        INTEGER(o_com_dfu)[h]  = community[h].day_last_followup;
        INTEGER(o_com_c2p)[h]  = community[h].c2p_group;
        INTEGER(o_com_eidx)[h] = (community[h].earliest_idx_day_ill == CB_MISSING) ?
                                  NA_INTEGER : community[h].earliest_idx_day_ill;
        INTEGER(o_com_ign)[h]  = community[h].ignore;
    }
    SEXP com_df, com_nm, com_rn;
    PROTECT(com_df = allocVector(VECSXP, 7)); n_protect++;
    PROTECT(com_nm = allocVector(STRSXP, 7)); n_protect++;
    SET_VECTOR_ELT(com_df, 0, o_com_id);   SET_STRING_ELT(com_nm, 0, mkChar("community"));
    SET_VECTOR_ELT(com_df, 1, o_com_ds);   SET_STRING_ELT(com_nm, 1, mkChar("day_epi_start"));
    SET_VECTOR_ELT(com_df, 2, o_com_dp);   SET_STRING_ELT(com_nm, 2, mkChar("day_epi_stop"));
    SET_VECTOR_ELT(com_df, 3, o_com_dfu);  SET_STRING_ELT(com_nm, 3, mkChar("day_last_followup"));
    SET_VECTOR_ELT(com_df, 4, o_com_c2p);  SET_STRING_ELT(com_nm, 4, mkChar("c2p_group"));
    SET_VECTOR_ELT(com_df, 5, o_com_eidx); SET_STRING_ELT(com_nm, 5, mkChar("earliest_idx_day_ill"));
    SET_VECTOR_ELT(com_df, 6, o_com_ign);  SET_STRING_ELT(com_nm, 6, mkChar("ignore"));
    setAttrib(com_df, R_NamesSymbol, com_nm);
    setAttrib(com_df, R_ClassSymbol, mkString("data.frame"));
    PROTECT(com_rn = allocVector(INTSXP, 2)); n_protect++;
    INTEGER(com_rn)[0] = NA_INTEGER; INTEGER(com_rn)[1] = -n_community;
    setAttrib(com_df, R_RowNamesSymbol, com_rn);

    /* time_dep_covariate: 4 columns, same number of rows as input */
    SEXP o_tdc_id, o_tdc_ds, o_tdc_dp, o_tdc_v, tdc_df, tdc_nm, tdc_rn;
    PROTECT(o_tdc_id = allocVector(INTSXP,  n_tdc_rows)); n_protect++;
    PROTECT(o_tdc_ds = allocVector(INTSXP,  n_tdc_rows)); n_protect++;
    PROTECT(o_tdc_dp = allocVector(INTSXP,  n_tdc_rows)); n_protect++;
    PROTECT(o_tdc_v  = allocVector(REALSXP, n_tdc_rows)); n_protect++;
    {
        /* Re-emit the same (id, day, day, updated_value) layout as the input.
           We reconstruct the layout from the updated people[] arrays. */
        int row = 0;
        for (i = 0; i < p_size; i++) {
            if (people[i].ignore != 0 || people[i].time_dep_covariate == NULL) {
                /* Person was ignored: emit the original input rows unchanged */
                /* Actually, ignored people still have rows in r_tdc from gen_population.
                   We copy them from r_tdc directly. */
                continue;
            }
            h = people[i].community;
            for (t = community[h].day_epi_start; t <= community[h].day_epi_stop; t++) {
                if (row >= n_tdc_rows) break;
                r = t - community[h].day_epi_start;
                INTEGER(o_tdc_id)[row] = i;
                INTEGER(o_tdc_ds)[row] = t;
                INTEGER(o_tdc_dp)[row] = t;
                REAL(o_tdc_v)[row]     = people[i].time_dep_covariate[r][0];
                row++;
            }
        }
        /* Fill any remaining rows with 0 (shouldn't happen for well-formed input) */
        for (; row < n_tdc_rows; row++) {
            INTEGER(o_tdc_id)[row] = NA_INTEGER;
            INTEGER(o_tdc_ds)[row] = NA_INTEGER;
            INTEGER(o_tdc_dp)[row] = NA_INTEGER;
            REAL(o_tdc_v)[row]     = 0.0;
        }
    }
    PROTECT(tdc_df = allocVector(VECSXP, 4)); n_protect++;
    PROTECT(tdc_nm = allocVector(STRSXP, 4)); n_protect++;
    SET_VECTOR_ELT(tdc_df, 0, o_tdc_id); SET_STRING_ELT(tdc_nm, 0, mkChar("id"));
    SET_VECTOR_ELT(tdc_df, 1, o_tdc_ds); SET_STRING_ELT(tdc_nm, 1, mkChar("day_start"));
    SET_VECTOR_ELT(tdc_df, 2, o_tdc_dp); SET_STRING_ELT(tdc_nm, 2, mkChar("day_stop"));
    SET_VECTOR_ELT(tdc_df, 3, o_tdc_v);  SET_STRING_ELT(tdc_nm, 3, mkChar("value"));
    setAttrib(tdc_df, R_NamesSymbol, tdc_nm);
    setAttrib(tdc_df, R_ClassSymbol, mkString("data.frame"));
    PROTECT(tdc_rn = allocVector(INTSXP, 2)); n_protect++;
    INTEGER(tdc_rn)[0] = NA_INTEGER; INTEGER(tdc_rn)[1] = -n_tdc_rows;
    setAttrib(tdc_df, R_RowNamesSymbol, tdc_rn);

    /* summary statistics */
    SEXP s_nidx, s_nsec, s_nsym, s_nasym, s_nesc, s_nimm;
    PROTECT(s_nidx  = ScalarInteger(n_sym_idx + n_asym_idx)); n_protect++;
    PROTECT(s_nsec  = ScalarInteger(n_sym_sec  + n_asym_sec)); n_protect++;
    PROTECT(s_nsym  = ScalarInteger(n_sym_sec)); n_protect++;
    PROTECT(s_nasym = ScalarInteger(n_asym_sec)); n_protect++;
    PROTECT(s_nesc  = ScalarInteger(n_esc_attacked)); n_protect++;
    PROTECT(s_nimm  = ScalarInteger(n_imm)); n_protect++;

    /* result list */
    SEXP result, res_nm;
    PROTECT(result = allocVector(VECSXP, 10)); n_protect++;
    PROTECT(res_nm = allocVector(STRSXP, 10)); n_protect++;
    SET_VECTOR_ELT(result, 0, pop_df);  SET_STRING_ELT(res_nm, 0, mkChar("pop"));
    SET_VECTOR_ELT(result, 1, tic_df);  SET_STRING_ELT(res_nm, 1, mkChar("time_ind_covariate"));
    SET_VECTOR_ELT(result, 2, com_df);  SET_STRING_ELT(res_nm, 2, mkChar("community"));
    SET_VECTOR_ELT(result, 3, tdc_df);  SET_STRING_ELT(res_nm, 3, mkChar("time_dep_covariate"));
    SET_VECTOR_ELT(result, 4, s_nidx);  SET_STRING_ELT(res_nm, 4, mkChar("n_index"));
    SET_VECTOR_ELT(result, 5, s_nsec);  SET_STRING_ELT(res_nm, 5, mkChar("n_secondary_inf"));
    SET_VECTOR_ELT(result, 6, s_nsym);  SET_STRING_ELT(res_nm, 6, mkChar("n_secondary_sym"));
    SET_VECTOR_ELT(result, 7, s_nasym); SET_STRING_ELT(res_nm, 7, mkChar("n_secondary_asym"));
    SET_VECTOR_ELT(result, 8, s_nesc);  SET_STRING_ELT(res_nm, 8, mkChar("n_escaped"));
    SET_VECTOR_ELT(result, 9, s_nimm);  SET_STRING_ELT(res_nm, 9, mkChar("n_preimmune"));
    setAttrib(result, R_NamesSymbol, res_nm);

    /* --- 14. Cleanup: free all C memory --- */

    free_arrays();  /* b, p, u, q, lb, lp, lu, lq, coeff_*, etc. */

    /* free community contact histories (linked lists) */
    for (h = 0; h < n_community; h++) {
        if (community[h].contact_history != NULL) {
            for (r = 0; r < community[h].epi_duration; r++) {
                ptr_contact = community[h].contact_history[r].c2p_contact;
                while (ptr_contact != NULL) {
                    ptr2_contact = ptr_contact->next;
                    free(ptr_contact); ptr_contact = ptr2_contact;
                }
                ptr_contact = community[h].contact_history[r].p2p_contact;
                while (ptr_contact != NULL) {
                    ptr2_contact = ptr_contact->next;
                    free(ptr_contact); ptr_contact = ptr2_contact;
                }
            }
            free(community[h].contact_history);
        }
        if (community[h].member != NULL) free(community[h].member);
        if (community[h].idx    != NULL) free(community[h].idx);
    }

    /* free per-person arrays */
    for (i = 0; i < p_size; i++) {
        if (people[i].time_ind_covariate != NULL) free(people[i].time_ind_covariate);
        if (people[i].time_dep_covariate != NULL) free_2d_array_double(people[i].time_dep_covariate);
        if (people[i].pat_covariate      != NULL) free_2d_array_double(people[i].pat_covariate);
        if (people[i].imm_covariate      != NULL) free(people[i].imm_covariate);
    }
    free(people);   people    = NULL;
    free(community); community = NULL;

    /* free cfg_pars dynamic arrays */
    if (cfg_pars.prob_incubation    != NULL) free(cfg_pars.prob_incubation);
    if (cfg_pars.prob_infectious    != NULL) free(cfg_pars.prob_infectious);
    if (cfg_pars.par_equiclass != NULL) {
        for (k = 0; k < cfg_pars.n_par_equiclass; k++)
            if (cfg_pars.par_equiclass[k].member != NULL) free(cfg_pars.par_equiclass[k].member);
        free(cfg_pars.par_equiclass);
    }
    if (cfg_pars.sim_par_effective  != NULL) free(cfg_pars.sim_par_effective);
    if (cfg_pars.par_fixed_id       != NULL) free(cfg_pars.par_fixed_id);
    if (cfg_pars.par_fixed_value    != NULL) free(cfg_pars.par_fixed_value);
    if (cfg_pars.c2p_covariate      != NULL) free(cfg_pars.c2p_covariate);
    if (cfg_pars.sus_p2p_covariate  != NULL) free(cfg_pars.sus_p2p_covariate);
    if (cfg_pars.inf_p2p_covariate  != NULL) free(cfg_pars.inf_p2p_covariate);
    if (cfg_pars.interaction        != NULL) free_2d_array_int(cfg_pars.interaction);
    if (cfg_pars.pat_covariate      != NULL) free(cfg_pars.pat_covariate);
    if (cfg_pars.imm_covariate      != NULL) free(cfg_pars.imm_covariate);

    UNPROTECT(n_protect);
    return result;
}

/* =========================================================
 * r_estimate_single
 *
 * Maximum-likelihood estimation for one epidemic data set
 * supplied as R data frames.  No files are read or written;
 * all population, contact, and imputation data come from R.
 *
 * Arguments:
 *   r_pop       – pop data frame (≥13 cols; cols 0-12 used)
 *   r_tic       – time_ind_covariate data frame (id + n_tic value cols),
 *                 or R_NilValue when n_time_ind_covariate == 0
 *   r_community – community data frame (5 cols)
 *   r_tdc       – time_dep_covariate data frame (id, day_start, day_stop,
 *                 + n_tdc value cols), or R_NilValue when n_tdc == 0
 *   r_c2p       – c2p_contact data frame or R_NilValue (auto-generate)
 *   r_p2p       – p2p_contact data frame or R_NilValue (auto-generate)
 *   r_impute    – impute data frame (9 cols) or R_NilValue (no EM)
 *   r_cfg       – named list from read_config() extended with transformed
 *                 fields (ini_par_effective_flat, lower_search_bound_tr,
 *                 upper_search_bound_tr, par_fixed_value_tr)
 *   r_i_inc     – 1-based incubation-period group index
 *   r_i_inf     – 1-based infectious-period group index
 *   r_seed      – integer RNG seed (used for EM state initialisation)
 *
 * Returns a named list:
 *   $est          – real vector, length n_par (raw scale: probs / ORs)
 *   $var_logit    – real vector, length n_par² (column-major covariance
 *                   matrix in logit/log scale; all-zeros if skip_variance)
 *   $log_likelihood – double scalar
 *   $error_code   – integer (0=OK, 1=max-eval, 2=var-error, -1=pre-fail)
 * ========================================================= */
SEXP r_estimate_single(SEXP r_pop, SEXP r_tic, SEXP r_community, SEXP r_tdc,
                       SEXP r_c2p, SEXP r_p2p, SEXP r_impute,
                       SEXP r_cfg, SEXP r_i_inc, SEXP r_i_inf, SEXP r_seed)
{
    /* ---- local variables ---- */
    int h, i, j, k, m, r, t;
    int i_inc, i_inf;
    int loc_p_size, loc_n_community;
    int len_inc, len_inf;
    int n_par, n_peq, n_ini;
    int n_b, n_p, n_u, n_q, n_c2p_cov, n_p2p_cov, n_pat_cov, n_imm_cov;
    int start_day, stop_day, contact_mode, ignore_flag;
    double offset_val, infective_prob;
    double log_likelihood;
    int estimation_error_type;
    int n_protect = 0;
    CONTACT *ptr_contact, *ptr2_contact, *ptr1_contact;
    CONTACT *sus_ptr_contact, *inf_ptr_contact;
    RISK    *ptr_risk, *ptr2_risk;
    RISK_CLASS *ptr_class, *ptr2_class;
    INTEGER_CHAIN *ptr_integer, *ptr2_integer;
    PEOPLE  *person, *member;
    STATE   *ptr_stat, *ptr2_stat;
    MATRIX   var, var_logit, der_mat;
    double  *est, *der, *der1, *der2, *value, *p2p_covariate_buf;

    /* ---- 1. Extract indices and seed ---- */
    i_inc = INTEGER(r_i_inc)[0] - 1;   /* convert to 0-based */
    i_inf = INTEGER(r_i_inf)[0] - 1;
    seed  = (long)INTEGER(r_seed)[0];   /* global; used later for EM */

    loc_p_size      = length(VECTOR_ELT(r_pop, 0));
    loc_n_community = length(VECTOR_ELT(r_community, 0));

    /* ---- 2. Build cfg_pars from r_cfg ---- */
    memset(&cfg_pars, 0, sizeof(CFG_PARS));

    /* basic counts */
    cfg_pars.n_b_mode             = cb_gi(r_cfg, "n_b_mode");
    cfg_pars.n_p_mode             = cb_gi(r_cfg, "n_p_mode");
    cfg_pars.n_u_mode             = cb_gi(r_cfg, "n_u_mode");
    cfg_pars.n_q_mode             = cb_gi(r_cfg, "n_q_mode");
    cfg_pars.n_time_ind_covariate = cb_gi(r_cfg, "n_time_ind_covariate");
    cfg_pars.n_time_dep_covariate = cb_gi(r_cfg, "n_time_dep_covariate");
    cfg_pars.n_covariate          = cfg_pars.n_time_ind_covariate + cfg_pars.n_time_dep_covariate;
    cfg_pars.n_c2p_covariate      = cb_gi(r_cfg, "n_c2p_covariate");
    cfg_pars.n_sus_p2p_covariate  = cb_gi(r_cfg, "n_sus_p2p_covariate");
    cfg_pars.n_inf_p2p_covariate  = cb_gi(r_cfg, "n_inf_p2p_covariate");
    cfg_pars.n_int_p2p_covariate  = cb_gi(r_cfg, "n_int_p2p_covariate");
    cfg_pars.n_p2p_covariate      = cb_gi(r_cfg, "n_p2p_covariate");
    cfg_pars.n_pat_covariate      = cb_gi(r_cfg, "n_pat_covariate");
    cfg_pars.n_imm_covariate      = cb_gi(r_cfg, "n_imm_covariate");
    cfg_pars.n_par                = cb_gi(r_cfg, "n_par");
    cfg_pars.n_par_equiclass      = cb_gi(r_cfg, "n_par_equiclass");
    cfg_pars.n_par_fixed          = cb_gi(r_cfg, "n_par_fixed");
    n_par     = cfg_pars.n_par;
    n_peq     = cfg_pars.n_par_equiclass;
    n_b       = cfg_pars.n_b_mode;
    n_p       = cfg_pars.n_p_mode;
    n_u       = cfg_pars.n_u_mode;
    n_q       = cfg_pars.n_q_mode;
    n_c2p_cov = cfg_pars.n_c2p_covariate;
    n_p2p_cov = cfg_pars.n_p2p_covariate;
    n_pat_cov = cfg_pars.n_pat_covariate;
    n_imm_cov = cfg_pars.n_imm_covariate;

    /* incubation period: select i_inc-th group */
    {
        SEXP min_inc_r  = cb_get(r_cfg, "min_inc");
        SEXP max_inc_r  = cb_get(r_cfg, "max_inc");
        SEXP prob_inc_r = cb_get(r_cfg, "prob_inc");
        cfg_pars.min_incubation = INTEGER(min_inc_r)[i_inc];
        cfg_pars.max_incubation = INTEGER(max_inc_r)[i_inc];
        len_inc = cfg_pars.max_incubation - cfg_pars.min_incubation + 1;
        cfg_pars.prob_incubation = (double *)malloc((size_t)(len_inc * sizeof(double)));
        SEXP pvec = VECTOR_ELT(prob_inc_r, i_inc);
        for (k = 0; k < len_inc; k++) cfg_pars.prob_incubation[k] = REAL(pvec)[k];
    }

    /* infectious period: select i_inf-th group */
    {
        SEXP lower_inf_r = cb_get(r_cfg, "lower_inf");
        SEXP upper_inf_r = cb_get(r_cfg, "upper_inf");
        SEXP prob_inf_r  = cb_get(r_cfg, "prob_inf");
        cfg_pars.lower_infectious = INTEGER(lower_inf_r)[i_inf];
        cfg_pars.upper_infectious = INTEGER(upper_inf_r)[i_inf];
        len_inf = cfg_pars.upper_infectious - cfg_pars.lower_infectious + 1;
        cfg_pars.prob_infectious = (double *)malloc((size_t)(len_inf * sizeof(double)));
        SEXP pvec = VECTOR_ELT(prob_inf_r, i_inf);
        for (k = 0; k < len_inf; k++) cfg_pars.prob_infectious[k] = REAL(pvec)[k];
        /* adjust for pre-onset infectiousness */
        if (cfg_pars.lower_infectious < 0) {
            int l2 = len_inc;
            double *pdf_inc = (double *)malloc((size_t)(l2 * sizeof(double)));
            double *sdf_inc = (double *)malloc((size_t)(l2 * sizeof(double)));
            for (k = 0; k < l2; k++) pdf_inc[k] = cfg_pars.prob_incubation[k];
            sdf_inc[0] = 1.0 - pdf_inc[0];
            for (k = 1; k < l2; k++)
                sdf_inc[k] = (k == l2-1) ? 0.0 : sdf_inc[k-1] - pdf_inc[k];
            for (t = cfg_pars.lower_infectious; t < 0; t++) {
                r = t - cfg_pars.lower_infectious;
                int linc = -t;
                if (linc >= cfg_pars.min_incubation) {
                    if (linc <= cfg_pars.max_incubation)
                        cfg_pars.prob_infectious[r] *= sdf_inc[linc - cfg_pars.min_incubation]
                                                     + pdf_inc[linc - cfg_pars.min_incubation];
                    else
                        cfg_pars.prob_infectious[r] = 0.0;
                }
            }
            free(pdf_inc); free(sdf_inc);
        }
    }

    /* covariate index arrays */
    {
        SEXP e; int n;
        n = cfg_pars.n_c2p_covariate;
        if (n > 0) {
            cfg_pars.c2p_covariate = (int *)malloc((size_t)(n * sizeof(int)));
            e = cb_get(r_cfg, "c2p_covariate");
            for (k = 0; k < n; k++) cfg_pars.c2p_covariate[k] = INTEGER(e)[k];
        }
        n = cfg_pars.n_sus_p2p_covariate;
        if (n > 0) {
            cfg_pars.sus_p2p_covariate = (int *)malloc((size_t)(n * sizeof(int)));
            e = cb_get(r_cfg, "sus_p2p_covariate");
            for (k = 0; k < n; k++) cfg_pars.sus_p2p_covariate[k] = INTEGER(e)[k];
        }
        n = cfg_pars.n_inf_p2p_covariate;
        if (n > 0) {
            cfg_pars.inf_p2p_covariate = (int *)malloc((size_t)(n * sizeof(int)));
            e = cb_get(r_cfg, "inf_p2p_covariate");
            for (k = 0; k < n; k++) cfg_pars.inf_p2p_covariate[k] = INTEGER(e)[k];
        }
        n = cfg_pars.n_int_p2p_covariate;
        if (n > 0) {
            make_2d_array_int(&cfg_pars.interaction, n, 2, 0);
            SEXP rint = cb_get(r_cfg, "interaction");
            for (k = 0; k < n; k++) {
                SEXP pair = VECTOR_ELT(rint, k);
                cfg_pars.interaction[k][0] = INTEGER(pair)[0];
                cfg_pars.interaction[k][1] = INTEGER(pair)[1];
            }
        }
        n = cfg_pars.n_pat_covariate;
        if (n > 0) {
            cfg_pars.pat_covariate = (int *)malloc((size_t)(n * sizeof(int)));
            e = cb_get(r_cfg, "pat_covariate");
            for (k = 0; k < n; k++) cfg_pars.pat_covariate[k] = INTEGER(e)[k];
        }
        n = cfg_pars.n_imm_covariate;
        if (n > 0) {
            cfg_pars.imm_covariate = (int *)malloc((size_t)(n * sizeof(int)));
            e = cb_get(r_cfg, "imm_covariate");
            for (k = 0; k < n; k++) cfg_pars.imm_covariate[k] = INTEGER(e)[k];
        }
    }

    /* parameter equiclasses */
    {
        cfg_pars.par_equiclass = NULL;
        if (n_peq > 0) {
            cfg_pars.par_equiclass = (PAR_EQUICLASS *)malloc((size_t)(n_peq * sizeof(PAR_EQUICLASS)));
            SEXP r_peq = cb_get(r_cfg, "par_equiclass");
            for (k = 0; k < n_peq; k++) {
                SEXP cls = VECTOR_ELT(r_peq, k);
                int sz   = cb_gi(cls, "size");
                cfg_pars.par_equiclass[k].size   = sz;
                cfg_pars.par_equiclass[k].member = (int *)malloc((size_t)(sz * sizeof(int)));
                SEXP memb = cb_get(cls, "member");
                for (j = 0; j < sz; j++)
                    cfg_pars.par_equiclass[k].member[j] = INTEGER(memb)[j];
            }
        }
    }

    /* sim_par_effective: not used for estimation; set NULL */
    cfg_pars.sim_par_effective = NULL;

    /* initial estimates (logit/log-transformed, row-major flat vector) */
    n_ini = cb_gi(r_cfg, "n_ini");
    cfg_pars.n_ini           = n_ini;
    cfg_pars.ini_par_provided = cb_gi(r_cfg, "ini_par_provided");
    cfg_pars.ini_par_effective = NULL;
    if (n_peq > 0) {
        make_2d_array_double(&cfg_pars.ini_par_effective, n_ini, n_peq, 0.0);
        SEXP ini_flat = cb_get(r_cfg, "ini_par_effective_flat");
        for (i = 0; i < n_ini; i++)
            for (j = 0; j < n_peq; j++)
                cfg_pars.ini_par_effective[i][j] = REAL(ini_flat)[i * n_peq + j];
    }

    /* search bounds (logit/log-transformed) */
    cfg_pars.search_bound_provided = cb_gi(r_cfg, "search_bound_provided");
    cfg_pars.lower_search_bound = cfg_pars.upper_search_bound = NULL;
    if (n_peq > 0) {
        make_1d_array_double(&cfg_pars.lower_search_bound, n_peq, 0.0);
        make_1d_array_double(&cfg_pars.upper_search_bound, n_peq, 0.0);
        SEXP lb_r = cb_get(r_cfg, "lower_search_bound_tr");
        SEXP ub_r = cb_get(r_cfg, "upper_search_bound_tr");
        for (k = 0; k < n_peq; k++) {
            cfg_pars.lower_search_bound[k] = REAL(lb_r)[k];
            cfg_pars.upper_search_bound[k] = REAL(ub_r)[k];
        }
    }

    /* convergence criteria (raw tolerances, not transformed) */
    cfg_pars.converge_criteria_provided = cb_gi(r_cfg, "converge_criteria_provided");
    cfg_pars.converge_criteria = NULL;
    if (n_peq > 0) {
        make_1d_array_double(&cfg_pars.converge_criteria, n_peq, 1e-8);
        if (cfg_pars.converge_criteria_provided == 1) {
            SEXP cc_r = cb_get(r_cfg, "converge_criteria");
            for (k = 0; k < n_peq; k++) cfg_pars.converge_criteria[k] = REAL(cc_r)[k];
        }
    }

    /* fixed parameters */
    cfg_pars.par_fixed_id    = NULL;
    cfg_pars.par_fixed_value = NULL;
    if (cfg_pars.n_par_fixed > 0) {
        cfg_pars.par_fixed_id    = (int *)malloc((size_t)(cfg_pars.n_par_fixed * sizeof(int)));
        cfg_pars.par_fixed_value = (double *)malloc((size_t)(cfg_pars.n_par_fixed * sizeof(double)));
        SEXP fid = cb_get(r_cfg, "par_fixed_id");
        SEXP fvl = cb_get(r_cfg, "par_fixed_value_tr");
        for (k = 0; k < cfg_pars.n_par_fixed; k++) {
            cfg_pars.par_fixed_id[k]    = INTEGER(fid)[k];
            cfg_pars.par_fixed_value[k] = REAL(fvl)[k];
        }
    }

    /* symptom-induced covariate flags */
    cfg_pars.PreIllness_as_covariate = cb_gi(r_cfg, "PreIllness_as_covariate");
    cfg_pars.PreIllness_covariate_id = cb_gi(r_cfg, "PreIllness_covariate_id");
    cfg_pars.Illness_as_covariate    = cb_gi(r_cfg, "Illness_as_covariate");
    cfg_pars.Illness_covariate_id    = cb_gi(r_cfg, "Illness_covariate_id");
    cfg_pars.RxIllness_as_covariate  = cb_gi(r_cfg, "RxIllness_as_covariate");
    cfg_pars.RxIllness_covariate_id  = cb_gi(r_cfg, "RxIllness_covariate_id");
    cfg_pars.RxIllness_prob          = cb_gd(r_cfg, "RxIllness_prob");
    cfg_pars.RxIllness_duration      = cb_gi(r_cfg, "RxIllness_duration");
    cfg_pars.RxIllness_index_only    = cb_gi(r_cfg, "RxIllness_index_only");

    /* contact history flags */
    cfg_pars.generate_c2p_contact                    = cb_gi(r_cfg, "generate_c2p_contact");
    cfg_pars.generate_p2p_contact                    = cb_gi(r_cfg, "generate_p2p_contact");
    cfg_pars.common_contact_history_within_community = cb_gi(r_cfg, "common_contact_history_within_community");
    cfg_pars.c2p_offset                             = cb_gi(r_cfg, "c2p_offset");
    cfg_pars.p2p_offset                             = cb_gi(r_cfg, "p2p_offset");

    /* estimation control flags */
    cfg_pars.simulation          = 0;
    cfg_pars.n_simulation        = 1;
    cfg_pars.simulation_only     = 0;
    cfg_pars.EM                  = cb_gi(r_cfg, "EM");
    cfg_pars.min_size_MCEM       = cb_gi(r_cfg, "min_size_MCEM");
    cfg_pars.n_base_sampling     = cb_gi(r_cfg, "n_base_sampling");
    cfg_pars.n_burnin_sampling   = cb_gi(r_cfg, "n_burnin_sampling");
    cfg_pars.n_burnin_iter       = cb_gi(r_cfg, "n_burnin_iter");
    cfg_pars.n_sampling_for_mce  = cb_gi(r_cfg, "n_sampling_for_mce");
    cfg_pars.use_bootstrap_for_mce = cb_gi(r_cfg, "use_bootstrap_for_mce");
    cfg_pars.skip_Evar_for_mce   = cb_gi(r_cfg, "skip_Evar_for_mce");
    cfg_pars.community_specific_weighting = cb_gi(r_cfg, "community_specific_weighting");
    cfg_pars.check_missingness   = 0;
    cfg_pars.check_mixing        = cb_gi(r_cfg, "check_mixing");
    cfg_pars.check_runtime       = cb_gi(r_cfg, "check_runtime");
    cfg_pars.optimization_choice = cb_gi(r_cfg, "optimization_choice");
    cfg_pars.preset_index        = cb_gi(r_cfg, "preset_index");
    cfg_pars.adjust_for_left_truncation  = cb_gi(r_cfg, "adjust_for_left_truncation");
    cfg_pars.adjust_for_right_censoring  = cb_gi(r_cfg, "adjust_for_right_censoring");
    cfg_pars.use_index_cases_to_improve_b = cb_gi(r_cfg, "use_index_cases_to_improve_b");
    cfg_pars.prop_mix_imm_esc    = cb_gd(r_cfg, "prop_mix_imm_esc");
    cfg_pars.asym_effect_est     = cb_gd(r_cfg, "asym_effect_est");
    cfg_pars.CPI_duration        = cb_gi(r_cfg, "CPI_duration");
    cfg_pars.skip_variance       = cb_gi(r_cfg, "skip_variance");
    cfg_pars.skip_output         = 1;   /* output handled by R, not files */
    cfg_pars.goodness_of_fit     = 0;
    cfg_pars.stat_test           = 0;
    cfg_pars.print_covariance    = 0;
    cfg_pars.silent_run          = cb_gi(r_cfg, "silent_run");
    cfg_pars.write_error_log     = 0;
    cfg_pars.idx_initiated_followup_in_simulation      = 0;
    cfg_pars.followup_duration_after_idx_in_simulation = 0;
    cfg_pars.SAR_n_covariate_sets       = 0;
    cfg_pars.SAR_sus_time_ind_covariate = NULL;
    cfg_pars.SAR_inf_time_ind_covariate = NULL;
    cfg_pars.SAR_sus_time_dep_covariate = NULL;
    cfg_pars.SAR_inf_time_dep_covariate = NULL;
    cfg_pars.n_inc               = cb_gi(r_cfg, "n_inc");
    cfg_pars.n_inf               = cb_gi(r_cfg, "n_inf");
    cfg_pars.n_c2p_group         = 0;
    cfg_pars.c2p_group           = NULL;

    /* effective infectious bounds (for SAR / R0 computation in R wrapper) */
    cfg_pars.effective_lower_infectious = NULL;
    cfg_pars.effective_upper_infectious = NULL;
    if (n_p > 0) {
        make_1d_array_int(&cfg_pars.effective_lower_infectious, n_p, 0);
        make_1d_array_int(&cfg_pars.effective_upper_infectious, n_p, 0);
        SEXP el_r = cb_get(r_cfg, "effective_lower_infectious_tr");
        SEXP eu_r = cb_get(r_cfg, "effective_upper_infectious_tr");
        for (k = 0; k < n_p; k++) {
            cfg_pars.effective_lower_infectious[k] = INTEGER(el_r)[k];
            cfg_pars.effective_upper_infectious[k] = INTEGER(eu_r)[k];
        }
    }

    /* R0 multiplier */
    cfg_pars.R0_multiplier_provided = cb_gi(r_cfg, "R0_multiplier_provided");
    cfg_pars.R0_multiplier          = NULL;
    cfg_pars.R0_multiplier_var      = NULL;
    if (cfg_pars.R0_multiplier_provided > 0 && n_p > 0) {
        make_1d_array_double(&cfg_pars.R0_multiplier,     n_p, 0.0);
        make_1d_array_double(&cfg_pars.R0_multiplier_var, n_p, 0.0);
        SEXP rm_r  = cb_get(r_cfg, "R0_multiplier");
        SEXP rmv_r = cb_get(r_cfg, "R0_multiplier_var");
        for (k = 0; k < n_p; k++) {
            cfg_pars.R0_multiplier[k]     = REAL(rm_r)[k];
            cfg_pars.R0_multiplier_var[k] = REAL(rmv_r)[k];
        }
    }

    /* ---- 3. Initialise globals ---- */
    p_size             = loc_p_size;
    n_community        = loc_n_community;
    importance_weight  = NULL;
    pdf_incubation = sdf_incubation = NULL;
    size_sample_states = 0;
    max_epi_duration   = 0;
    ITER               = 0;
    n_need_OEM         = 0;
    n_need_MCEM        = 0;

    /* ---- 4. Initialise MATRIX objects ---- */
    initialize_matrix(&var);
    initialize_matrix(&var_logit);
    inflate_matrix(&var,      n_par, n_par, 0.0);
    inflate_matrix(&var_logit, n_par, n_par, 0.0);
    initialize_matrix(&der_mat);
    inflate_matrix(&der_mat, n_p, n_par, 0.0);

    /* ---- 5. pdf_incubation / sdf_incubation (core.h lines 70-107) ---- */
    {
        int l2 = len_inc;
        pdf_incubation = (double *)malloc((size_t)(l2 * sizeof(double)));
        sdf_incubation = (double *)malloc((size_t)(l2 * sizeof(double)));
        for (t = cfg_pars.min_incubation; t <= cfg_pars.max_incubation; t++) {
            r = t - cfg_pars.min_incubation;
            pdf_incubation[r] = cfg_pars.prob_incubation[r];
        }
        for (t = cfg_pars.min_incubation; t <= cfg_pars.max_incubation; t++) {
            r = t - cfg_pars.min_incubation;
            if (t == cfg_pars.min_incubation)       sdf_incubation[r] = 1.0 - pdf_incubation[r];
            else if (t == cfg_pars.max_incubation)  sdf_incubation[r] = 0.0;
            else                                     sdf_incubation[r] = sdf_incubation[r-1] - pdf_incubation[r];
        }
        /* apply sdf to pre-onset infectious probability */
        for (t = cfg_pars.lower_infectious; t <= cfg_pars.upper_infectious; t++) {
            r = t - cfg_pars.lower_infectious;
            if (t < 0) {
                int linc = -t;
                if (linc >= cfg_pars.min_incubation) {
                    if (linc <= cfg_pars.max_incubation)
                        cfg_pars.prob_infectious[r] *=
                            sdf_incubation[linc - cfg_pars.min_incubation] +
                            pdf_incubation[linc - cfg_pars.min_incubation];
                    else
                        cfg_pars.prob_infectious[r] = 0.0;
                }
            }
        }
    }

    /* ---- 6. Working arrays ---- */
    est = der = der1 = der2 = value = p2p_covariate_buf = NULL;
    if (n_par > 0) {
        est  = (double *)malloc((size_t)(n_par * sizeof(double)));
        der  = (double *)malloc((size_t)(n_par * sizeof(double)));
        der1 = (double *)malloc((size_t)(n_par * sizeof(double)));
        der2 = (double *)malloc((size_t)(n_par * sizeof(double)));
    }
    {
        int k2 = cfg_pars.n_time_ind_covariate > cfg_pars.n_time_dep_covariate ?
                 cfg_pars.n_time_ind_covariate : cfg_pars.n_time_dep_covariate;
        if (k2 > 0) value = (double *)malloc((size_t)(k2 * sizeof(double)));
    }
    if (n_p2p_cov > 0)
        p2p_covariate_buf = (double *)malloc((size_t)(n_p2p_cov * sizeof(double)));

    /* ---- 7. Allocate and fill people[] ---- */
    people = (PEOPLE *)malloc((size_t)(p_size * sizeof(PEOPLE)));
    if (!people) error("r_estimate_single: cannot allocate people array");
    {
        int    *p_id   = INTEGER(VECTOR_ELT(r_pop,  0));
        int    *p_comm = INTEGER(VECTOR_ELT(r_pop,  1));
        int    *p_pimm = INTEGER(VECTOR_ELT(r_pop,  2));
        int    *p_inf  = INTEGER(VECTOR_ELT(r_pop,  3));
        int    *p_sym  = INTEGER(VECTOR_ELT(r_pop,  4));
        int    *p_dill = INTEGER(VECTOR_ELT(r_pop,  5));
        int    *p_exit = INTEGER(VECTOR_ELT(r_pop,  6));
        int    *p_dex  = INTEGER(VECTOR_ELT(r_pop,  7));
        int    *p_idx  = INTEGER(VECTOR_ELT(r_pop,  8));
        int    *p_umod = INTEGER(VECTOR_ELT(r_pop,  9));
        int    *p_qmod = INTEGER(VECTOR_ELT(r_pop, 10));
        double *p_wt   = REAL   (VECTOR_ELT(r_pop, 11));
        int    *p_ign  = INTEGER(VECTOR_ELT(r_pop, 12));

        for (i = 0; i < p_size; i++) {
            people[i].id         = p_id[i];
            people[i].community  = p_comm[i];
            people[i].pre_immune = p_pimm[i];
            people[i].infection  = p_inf[i];
            people[i].symptom    = p_sym[i];
            people[i].day_ill    = p_dill[i];
            people[i].exit       = p_exit[i];
            people[i].day_exit   = p_dex[i];
            people[i].idx        = p_idx[i];
            people[i].u_mode     = p_umod[i];
            people[i].q_mode     = p_qmod[i];
            people[i].weight     = p_wt[i];
            people[i].ignore     = p_ign[i];
            people[i].day_infection       = CB_MISSING;
            people[i].day_infection_lower = CB_MISSING;
            people[i].day_infection_upper = CB_MISSING;
            people[i].day_infective_lower = CB_MISSING;
            people[i].day_infective_upper = CB_MISSING;
            people[i].final_risk_day      = CB_MISSING;
            people[i].any_exposure        = 0;
            people[i].current_state       = 0;
            people[i].size_possible_states = 0;
            people[i].time_ind_covariate  = NULL;
            people[i].time_dep_covariate  = NULL;
            people[i].contact_history     = NULL;
            people[i].risk_history        = NULL;
            people[i].risk_class          = NULL;
            people[i].imm_covariate       = NULL;
            people[i].pat_covariate       = NULL;
            people[i].possible_states     = NULL;
            /* fix sentinels for EM==0: no infection → day_ill = MISSING */
            if (people[i].infection == 0) people[i].day_ill = CB_MISSING;
            if (people[i].exit == 0)      people[i].day_exit = 1000000;
        }
    }

    /* ---- 8. Allocate community[] ---- */
    community = (COMMUNITY *)malloc((size_t)(n_community * sizeof(COMMUNITY)));
    if (!community) { free(people); error("r_estimate_single: cannot allocate community array"); }
    {
        int *c_id  = INTEGER(VECTOR_ELT(r_community, 0));
        int *c_ds  = INTEGER(VECTOR_ELT(r_community, 1));
        int *c_dp  = INTEGER(VECTOR_ELT(r_community, 2));
        int *c_dfu = INTEGER(VECTOR_ELT(r_community, 3));
        int *c_c2p = INTEGER(VECTOR_ELT(r_community, 4));

        for (h = 0; h < n_community; h++) {
            community[h].id                   = c_id[h];
            community[h].day_epi_start        = c_ds[h];
            community[h].day_epi_stop         = c_dp[h];
            community[h].day_last_followup    = c_dfu[h];
            community[h].c2p_group            = c_c2p[h];
            community[h].epi_duration         = c_dp[h] - c_ds[h] + 1;
            if (community[h].epi_duration > max_epi_duration)
                max_epi_duration = community[h].epi_duration;
            community[h].size                 = 0;
            community[h].size_idx             = 0;
            community[h].earliest_idx_day_ill = CB_MISSING;
            community[h].latest_idx_day_ill   = CB_MISSING;
            community[h].counter              = 0;
            community[h].ignore               = 0;
            community[h].member               = NULL;
            community[h].idx                  = NULL;
            community[h].member_impute        = NULL;
            community[h].size_impute          = 0;
            community[h].size_possible_states = 1;
            community[h].contact_history      = NULL;
            community[h].risk_class           = NULL;
            community[h].risk_class_rear      = NULL;
            community[h].sample_states        = NULL;
            community[h].sample_states_rear   = NULL;
            community[h].list_states          = NULL;
            community[h].list_states_rear     = NULL;
        }
    }

    /* build community membership */
    for (i = 0; i < p_size; i++)
        if (people[i].ignore == 0) community[people[i].community].size++;
    for (h = 0; h < n_community; h++)
        if (community[h].size > 0)
            community[h].member = (int *)malloc((size_t)(community[h].size * sizeof(int)));
    for (h = 0; h < n_community; h++) community[h].counter = 0;
    for (i = 0; i < p_size; i++)
        if (people[i].ignore == 0) {
            h = people[i].community;
            community[h].member[community[h].counter++] = i;
        }

    /* ---- 9. Allocate per-person covariate arrays ---- */
    {
        int n_tic = cfg_pars.n_time_ind_covariate;
        int n_tdc = cfg_pars.n_time_dep_covariate;
        for (i = 0; i < p_size; i++) {
            if (people[i].ignore == 0) {
                h = people[i].community;
                if (n_tic > 0)
                    make_1d_array_double(&people[i].time_ind_covariate, n_tic, 0.0);
                if (n_tdc > 0)
                    make_2d_array_double(&people[i].time_dep_covariate,
                                         community[h].epi_duration, n_tdc, 0.0);
            }
        }
    }

    /* fill time-independent covariates (all n_tic value columns) */
    if (cfg_pars.n_time_ind_covariate > 0 && r_tic != R_NilValue) {
        int n_tic_rows = length(VECTOR_ELT(r_tic, 0));
        int *tic_id    = INTEGER(VECTOR_ELT(r_tic, 0));
        for (m = 0; m < n_tic_rows; m++) {
            i = tic_id[m];
            if (i >= 0 && i < p_size && people[i].ignore == 0 &&
                people[i].time_ind_covariate != NULL) {
                for (j = 0; j < cfg_pars.n_time_ind_covariate; j++) {
                    double *col = REAL(VECTOR_ELT(r_tic, j + 1));
                    people[i].time_ind_covariate[j] = col[m];
                }
            }
        }
    }

    /* fill time-dependent covariates (all n_tdc value columns) */
    if (cfg_pars.n_time_dep_covariate > 0 && r_tdc != R_NilValue) {
        int n_tdc_rows = length(VECTOR_ELT(r_tdc, 0));
        int *tdc_id    = INTEGER(VECTOR_ELT(r_tdc, 0));
        int *tdc_ds    = INTEGER(VECTOR_ELT(r_tdc, 1));
        int *tdc_dp    = INTEGER(VECTOR_ELT(r_tdc, 2));
        for (m = 0; m < n_tdc_rows; m++) {
            i = tdc_id[m];
            if (i < 0 || i >= p_size || people[i].ignore != 0 ||
                people[i].time_dep_covariate == NULL) continue;
            h = people[i].community;
            for (t = tdc_ds[m]; t <= tdc_dp[m]; t++) {
                if (t >= community[h].day_epi_start && t <= community[h].day_epi_stop) {
                    r = t - community[h].day_epi_start;
                    for (j = 0; j < cfg_pars.n_time_dep_covariate; j++) {
                        double *col = REAL(VECTOR_ELT(r_tdc, j + 3));
                        people[i].time_dep_covariate[r][j] = col[m];
                    }
                }
            }
        }
    }

    /* ---- 9b. Infection-day ranges for EM==0 (mirrors core.h lines 388-404) ---- */
    if (cfg_pars.EM == 0) {
        for (i = 0; i < p_size; i++) {
            if (people[i].ignore == 0 && people[i].infection == 1) {
                h = people[i].community;
                if (people[i].day_ill > community[h].day_epi_stop &&
                    cfg_pars.adjust_for_right_censoring == 1) {
                    people[i].infection = 0;
                    people[i].symptom   = 0;
                    people[i].day_ill   = CB_MISSING;
                    continue;
                }
                people[i].day_infection_lower =
                    CB_max(people[i].day_ill - cfg_pars.max_incubation, community[h].day_epi_start);
                people[i].day_infection_upper =
                    CB_max(people[i].day_ill - cfg_pars.min_incubation, community[h].day_epi_start);
                people[i].day_infection_upper =
                    CB_min(people[i].day_infection_upper, community[h].day_epi_stop);
                people[i].day_infective_lower = people[i].day_ill + cfg_pars.lower_infectious;
                people[i].day_infective_upper = people[i].day_ill + cfg_pars.upper_infectious;
            }
        }
    }

    /* ---- 10. Allocate contact history structures ---- */
    /* community-level (shared contact history) */
    if (cfg_pars.common_contact_history_within_community == 1) {
        for (h = 0; h < n_community; h++) {
            if (community[h].size > 0) {
                community[h].contact_history = (CONTACT_HISTORY *)malloc(
                    (size_t)(community[h].epi_duration * sizeof(CONTACT_HISTORY)));
                for (t = community[h].day_epi_start; t <= community[h].day_epi_stop; t++) {
                    r = t - community[h].day_epi_start;
                    community[h].contact_history[r].c2p_contact      = NULL;
                    community[h].contact_history[r].c2p_contact_rear = NULL;
                    community[h].contact_history[r].p2p_contact      = NULL;
                    community[h].contact_history[r].p2p_contact_rear = NULL;
                }
            }
        }
    } else {
        /* person-level contact history and risk history */
        for (i = 0; i < p_size; i++) {
            if (people[i].ignore == 0) {
                h = people[i].community;
                people[i].contact_history = (CONTACT_HISTORY *)malloc(
                    (size_t)(community[h].epi_duration * sizeof(CONTACT_HISTORY)));
                people[i].risk_history = (RISK_HISTORY *)malloc(
                    (size_t)(community[h].epi_duration * sizeof(RISK_HISTORY)));
                for (t = community[h].day_epi_start; t <= community[h].day_epi_stop; t++) {
                    r = t - community[h].day_epi_start;
                    people[i].contact_history[r].c2p_contact      = NULL;
                    people[i].contact_history[r].c2p_contact_rear = NULL;
                    people[i].contact_history[r].p2p_contact      = NULL;
                    people[i].contact_history[r].p2p_contact_rear = NULL;
                    people[i].risk_history[r].c2p_risk      = NULL;
                    people[i].risk_history[r].c2p_risk_rear = NULL;
                    people[i].risk_history[r].p2p_risk      = NULL;
                    people[i].risk_history[r].p2p_risk_rear = NULL;
                }
            }
        }
    }

    /* ---- 11. Global parameter arrays ---- */
    create_arrays();

    /* ---- 12. Build c2p contact history ---- */
    if (r_c2p != R_NilValue) {
        /* consistency check */
        if (cfg_pars.generate_c2p_contact == 1 && cfg_pars.silent_run == 0)
            Rprintf("Warning: c2p_contact provided but generate_c2p_contact=1; using provided data.\n");
        {
            int n_c2p_rows = length(VECTOR_ELT(r_c2p, 0));
            int *c2p_col0  = INTEGER(VECTOR_ELT(r_c2p, 0)); /* comm_id or person_id */
            int *c2p_ds    = INTEGER(VECTOR_ELT(r_c2p, 1));
            int *c2p_dp    = INTEGER(VECTOR_ELT(r_c2p, 2));
            int *c2p_cm    = INTEGER(VECTOR_ELT(r_c2p, 3));
            double *c2p_off = REAL(VECTOR_ELT(r_c2p, 4));
            int *c2p_ign   = INTEGER(VECTOR_ELT(r_c2p, 5));

            if (cfg_pars.common_contact_history_within_community == 1) {
                for (m = 0; m < n_c2p_rows; m++) {
                    h = c2p_col0[m];
                    offset_val   = (cfg_pars.c2p_offset == 0) ? 0.0 : c2p_off[m];
                    ignore_flag  = c2p_ign[m];
                    if (ignore_flag != 1 && h >= 0 && h < n_community && community[h].size > 0)
                        add_c2p_contact_history_to_community(h, c2p_ds[m], c2p_dp[m],
                                                              c2p_cm[m], offset_val);
                }
            } else {
                for (m = 0; m < n_c2p_rows; m++) {
                    i = c2p_col0[m];
                    offset_val  = (cfg_pars.c2p_offset == 0) ? 0.0 : c2p_off[m];
                    ignore_flag = c2p_ign[m];
                    if (ignore_flag != 1 && i >= 0 && i < p_size && people[i].ignore != 1) {
                        h = people[i].community;
                        for (t = c2p_ds[m]; t <= c2p_dp[m]; t++) {
                            if (t >= community[h].day_epi_start &&
                                t <= community[h].day_epi_stop  &&
                                t <= people[i].day_exit)
                                add_c2p_contact_history(t, people + i, c2p_cm[m], offset_val);
                        }
                    }
                }
            }
        }
    } else {
        /* auto-generate: random mixing */
        if (cfg_pars.generate_c2p_contact == 0 && cfg_pars.silent_run == 0)
            Rprintf("Warning: c2p_contact absent but generate_c2p_contact=0; auto-generating.\n");
        if (cfg_pars.common_contact_history_within_community == 1) {
            for (h = 0; h < n_community; h++)
                if (community[h].size > 0)
                    add_c2p_contact_history_to_community(h,
                        community[h].day_epi_start, community[h].day_epi_stop, 0, 0.0);
        } else {
            for (i = 0; i < p_size; i++)
                if (people[i].ignore == 0) {
                    h = people[i].community;
                    for (t = community[h].day_epi_start; t <= community[h].day_epi_stop; t++)
                        if (t <= people[i].day_exit)
                            add_c2p_contact_history(t, people + i, 0, 0.0);
                }
        }
    }

    /* ---- 13. Build p2p contact history ---- */
    if (r_p2p != R_NilValue) {
        if (cfg_pars.generate_p2p_contact == 1 && cfg_pars.silent_run == 0)
            Rprintf("Warning: p2p_contact provided but generate_p2p_contact=1; using provided data.\n");
        {
            int n_p2p_rows = length(VECTOR_ELT(r_p2p, 0));
            if (cfg_pars.common_contact_history_within_community == 1) {
                /* shared: community_id, start_day, stop_day, contact_mode, offset, ignore */
                int *p2p_col0  = INTEGER(VECTOR_ELT(r_p2p, 0));
                int *p2p_ds    = INTEGER(VECTOR_ELT(r_p2p, 1));
                int *p2p_dp    = INTEGER(VECTOR_ELT(r_p2p, 2));
                int *p2p_cm    = INTEGER(VECTOR_ELT(r_p2p, 3));
                double *p2p_off = REAL(VECTOR_ELT(r_p2p, 4));
                int *p2p_ign   = INTEGER(VECTOR_ELT(r_p2p, 5));
                for (m = 0; m < n_p2p_rows; m++) {
                    h = p2p_col0[m];
                    offset_val  = (cfg_pars.p2p_offset == 0) ? 0.0 : p2p_off[m];
                    if (p2p_ign[m] != 1 && h >= 0 && h < n_community && community[h].size > 0)
                        add_p2p_contact_history_to_community(h, p2p_ds[m], p2p_dp[m],
                                                              p2p_cm[m], offset_val);
                }
            } else {
                /* custom: start_day, stop_day, person_i, person_j, contact_mode, offset, ignore */
                int *p2p_ds    = INTEGER(VECTOR_ELT(r_p2p, 0));
                int *p2p_dp    = INTEGER(VECTOR_ELT(r_p2p, 1));
                int *p2p_pi    = INTEGER(VECTOR_ELT(r_p2p, 2));
                int *p2p_pj    = INTEGER(VECTOR_ELT(r_p2p, 3));
                int *p2p_cm    = INTEGER(VECTOR_ELT(r_p2p, 4));
                double *p2p_off = REAL(VECTOR_ELT(r_p2p, 5));
                int *p2p_ign   = INTEGER(VECTOR_ELT(r_p2p, 6));
                for (m = 0; m < n_p2p_rows; m++) {
                    i = p2p_pi[m]; j = p2p_pj[m];
                    offset_val = (cfg_pars.p2p_offset == 0) ? 0.0 : p2p_off[m];
                    if (p2p_ign[m] != 1 &&
                        i >= 0 && i < p_size && j >= 0 && j < p_size &&
                        people[i].ignore == 0 && people[j].ignore == 0 &&
                        people[i].community == people[j].community &&
                        people[i].id != people[j].id) {
                        h = people[i].community;
                        for (t = p2p_ds[m]; t <= p2p_dp[m]; t++) {
                            if (t >= community[h].day_epi_start &&
                                t <= community[h].day_epi_stop  &&
                                t <= people[i].day_exit && t <= people[j].day_exit) {
                                ptr1_contact = add_p2p_contact_history(
                                    t, people+i, people+j, p2p_cm[m], offset_val);
                                sus_ptr_contact = add_p2p_contact_history(
                                    t, people+j, people+i, p2p_cm[m], offset_val);
                                ptr1_contact->pair  = sus_ptr_contact;
                                sus_ptr_contact->pair = ptr1_contact;
                            }
                        }
                    }
                }
            }
        }
    } else {
        /* auto-generate */
        if (cfg_pars.generate_p2p_contact == 0 && cfg_pars.silent_run == 0)
            Rprintf("Warning: p2p_contact absent but generate_p2p_contact=0; auto-generating.\n");
        if (cfg_pars.common_contact_history_within_community == 1) {
            for (h = 0; h < n_community; h++)
                if (community[h].size > 0)
                    add_p2p_contact_history_to_community(h,
                        community[h].day_epi_start, community[h].day_epi_stop, 0, 0.0);
        }
        /* individualized p2p auto-generation is not supported */
    }

    /* ---- 14. Risk classes ---- */
    create_risk_class();

    /* ---- 15. Case-ascertained design: set index cases ---- */
    /* Always performed when adjust_for_left_truncation==1 (not conditional on
       simulation==0 as in core.h), because simulated data has already lost its
       C state after simulate_single() returned to R.                         */
    if (cfg_pars.adjust_for_left_truncation == 1) {
        for (h = 0; h < n_community; h++) {
            set_index_cases(community + h);
            if (!(cfg_pars.preset_index == 0 && cfg_pars.EM == 1)) {
                if (community[h].size_idx == 0 ||
                    (community[h].size_idx == community[h].size &&
                     cfg_pars.use_index_cases_to_improve_b == 0)) {
                    community[h].ignore = 1;
                    if (community[h].idx != NULL) {
                        free(community[h].idx);
                        community[h].idx = NULL;
                    }
                    for (j = 0; j < community[h].size; j++)
                        people[community[h].member[j]].ignore = 1;
                }
            }
        }
    }

    /* ---- 16. EM: build possible-states from impute data frame ---- */
    if (cfg_pars.EM == 1) {
        if (r_impute != R_NilValue) {
            int n_imp = length(VECTOR_ELT(r_impute, 0));
            int *imp_pid  = INTEGER(VECTOR_ELT(r_impute, 0));
            int *imp_pimm = INTEGER(VECTOR_ELT(r_impute, 1));
            int *imp_pesc = INTEGER(VECTOR_ELT(r_impute, 2));
            int *imp_psym = INTEGER(VECTOR_ELT(r_impute, 3));
            int *imp_pss  = INTEGER(VECTOR_ELT(r_impute, 4));
            int *imp_psp  = INTEGER(VECTOR_ELT(r_impute, 5));
            int *imp_pasy = INTEGER(VECTOR_ELT(r_impute, 6));
            int *imp_pass = INTEGER(VECTOR_ELT(r_impute, 7));
            int *imp_pasp = INTEGER(VECTOR_ELT(r_impute, 8));

            for (m = 0; m < n_imp; m++) {
                i = imp_pid[m];
                if (i < 0 || i >= p_size || people[i].ignore == 1) continue;
                {
                    int pimm = imp_pimm[m], pesc = imp_pesc[m];
                    int psym = imp_psym[m], pss  = imp_pss[m],  psp  = imp_psp[m];
                    int pasy = imp_pasy[m], pass = imp_pass[m], pasp = imp_pasp[m];
                    int sz = 0, kk;
                    if (pimm == 1) sz++;
                    if (pesc == 1) sz++;
                    if (psym == 1) sz += (psp - pss + 1);
                    if (pasy == 1) sz += (pasp - pass + 1);
                    people[i].size_possible_states = sz;
                    if (sz > 0) {
                        make_2d_array_int(&people[i].possible_states, 2, sz, 0);
                        kk = 0;
                        if (pimm == 1) { people[i].possible_states[0][kk] = -1000000; kk++; }
                        if (pesc == 1) { people[i].possible_states[0][kk] =  1000000; kk++; }
                        if (psym == 1) {
                            for (j = pss; j <= psp; j++) {
                                people[i].possible_states[0][kk] = j;
                                people[i].possible_states[1][kk] = 1;
                                kk++;
                            }
                        }
                        if (pasy == 1) {
                            for (j = pass; j <= pasp; j++) {
                                people[i].possible_states[0][kk] = j;
                                people[i].possible_states[1][kk] = 0;
                                kk++;
                            }
                        }
                    }
                }
            }
        }

        /* ---- 17. Count uncertain individuals per community ---- */
        for (h = 0; h < n_community; h++) {
            for (j = 0; j < community[h].size; j++) {
                i = community[h].member[j];
                if (people[i].size_possible_states > 0) {
                    community[h].size_impute++;
                    community[h].size_possible_states *= people[i].size_possible_states;
                }
            }
            if (community[h].size_impute > 0) {
                community[h].member_impute =
                    (int *)malloc((size_t)(community[h].size_impute * sizeof(int)));
                k = 0;
                for (j = 0; j < community[h].size; j++) {
                    i = community[h].member[j];
                    if (people[i].size_possible_states > 0)
                        community[h].member_impute[k++] = i;
                }
            }
            if (community[h].size_possible_states > 1 &&
                community[h].size_possible_states < cfg_pars.min_size_MCEM)
                n_need_OEM++;
            if (community[h].size_possible_states >= cfg_pars.min_size_MCEM)
                n_need_MCEM++;
        }

        /* ---- 18. Infection-day ranges for certain individuals ---- */
        for (i = 0; i < p_size; i++) {
            if (people[i].ignore == 0 && people[i].infection == 1 &&
                people[i].size_possible_states == 0) {
                h = people[i].community;
                if (people[i].day_ill > community[h].day_epi_stop &&
                    cfg_pars.adjust_for_right_censoring == 1) {
                    people[i].infection = 0;
                    people[i].symptom   = 0;
                    people[i].day_ill   = CB_MISSING;
                    continue;
                }
                people[i].day_infection_lower =
                    CB_max(people[i].day_ill - cfg_pars.max_incubation, community[h].day_epi_start);
                people[i].day_infection_upper =
                    CB_max(people[i].day_ill - cfg_pars.min_incubation, community[h].day_epi_start);
                people[i].day_infection_upper =
                    CB_min(people[i].day_infection_upper, community[h].day_epi_stop);
                people[i].day_infective_lower = people[i].day_ill + cfg_pars.lower_infectious;
                people[i].day_infective_upper = people[i].day_ill + cfg_pars.upper_infectious;
            }
        }

        /* ---- 19. Initialise uncertain states (draws from RNG) ---- */
        mt_init(seed);
        for (h = 0; h < n_community; h++) {
            if (community[h].size > 0) {
                for (m = 0; m < community[h].size; m++) {
                    person = people + community[h].member[m];
                    if (person->size_possible_states > 0) {
                        person->current_state =
                            (int)floor(person->size_possible_states * runiform(&seed));
                        k = person->current_state;
                        if (person->possible_states[0][k] == -1000000) {
                            person->pre_immune = 1;
                            person->infection  = 0;
                            person->symptom    = 0;
                            person->day_ill    = CB_MISSING;
                        } else if (person->possible_states[0][k] == 1000000) {
                            person->pre_immune = 0;
                            person->infection  = 0;
                            person->symptom    = 0;
                            person->day_ill    = CB_MISSING;
                        } else {
                            person->pre_immune = 0;
                            person->infection  = 1;
                            person->symptom    = person->possible_states[1][k];
                            person->day_ill    = person->possible_states[0][k];
                            person->day_infective_lower =
                                person->day_ill + cfg_pars.lower_infectious;
                            person->day_infective_upper =
                                person->day_ill + cfg_pars.upper_infectious;
                            person->day_infection_lower =
                                CB_max(person->day_ill - cfg_pars.max_incubation,
                                       community[h].day_epi_start);
                            person->day_infection_upper =
                                CB_max(person->day_ill - cfg_pars.min_incubation,
                                       community[h].day_epi_start);
                            person->day_infection_upper =
                                CB_min(person->day_infection_upper, community[h].day_epi_stop);
                        }
                    } else {
                        if (person->pre_immune == 1)
                            person->final_risk_day = -1000000;
                        else {
                            if (person->infection == 1)
                                person->final_risk_day =
                                    person->day_ill - cfg_pars.min_incubation;
                            else
                                person->final_risk_day = community[h].day_epi_stop;
                        }
                    }
                }
            }
        }

        /* re-set index cases after sampling uncertain states */
        if (cfg_pars.adjust_for_left_truncation == 1 && cfg_pars.preset_index == 0) {
            for (h = 0; h < n_community; h++)
                set_index_cases(community + h);
        }
    } /* end EM==1 block */

    /* ---- 20. Build risk history ---- */
    if (cfg_pars.common_contact_history_within_community == 1) {
        /* community-level risk history via risk classes */
        for (h = 0; h < n_community; h++) {
            if (community[h].ignore == 0 && community[h].size > 0) {
                for (t = community[h].day_epi_start; t <= community[h].day_epi_stop; t++) {
                    r = t - community[h].day_epi_start;
                    ptr_contact = community[h].contact_history[r].c2p_contact;
                    while (ptr_contact != NULL) {
                        add_c2p_risk_history_to_risk_class(
                            t, h, ptr_contact->contact_mode, ptr_contact->offset);
                        ptr_contact = ptr_contact->next;
                    }
                }
            }
        }
        for (h = 0; h < n_community; h++) {
            if (community[h].ignore == 0 && community[h].size > 0) {
                for (t = community[h].day_epi_start; t <= community[h].day_epi_stop; t++) {
                    r = t - community[h].day_epi_start;
                    ptr_contact = community[h].contact_history[r].p2p_contact;
                    while (ptr_contact != NULL) {
                        member = people + ptr_contact->contact_id;
                        if (member->infection == 1 &&
                            member->day_infective_lower != CB_MISSING &&
                            member->day_infective_upper != CB_MISSING &&
                            t >= member->day_infective_lower &&
                            t <= member->day_infective_upper &&
                            t <= member->day_exit) {
                            infective_prob =
                                cfg_pars.prob_infectious[t - member->day_infective_lower];
                            add_p2p_risk_history_to_risk_class(
                                t, member, ptr_contact->contact_mode,
                                ptr_contact->offset, infective_prob, member->symptom);
                        }
                        ptr_contact = ptr_contact->next;
                    }
                }
            }
        }
    } else {
        /* individual-level risk history */
        for (i = 0; i < p_size; i++) {
            person = people + i;
            if (person->ignore == 0) {
                h = person->community;
                /* c2p risk */
                for (t = community[h].day_epi_start; t <= community[h].day_epi_stop; t++) {
                    r = t - community[h].day_epi_start;
                    ptr_contact = person->contact_history[r].c2p_contact;
                    while (ptr_contact != NULL) {
                        add_c2p_risk_history(t, person,
                            ptr_contact->contact_mode, ptr_contact->offset);
                        ptr_contact = ptr_contact->next;
                    }
                }
                /* p2p risk */
                if (person->infection == 1 &&
                    person->day_infective_lower != CB_MISSING &&
                    person->day_infective_upper != CB_MISSING) {
                    for (t = person->day_infective_lower; t <= person->day_infective_upper; t++) {
                        if (t >= community[h].day_epi_start && t <= community[h].day_epi_stop) {
                            r = t - community[h].day_epi_start;
                            infective_prob =
                                cfg_pars.prob_infectious[t - person->day_infective_lower];
                            inf_ptr_contact = person->contact_history[r].p2p_contact;
                            while (inf_ptr_contact != NULL) {
                                member = people + inf_ptr_contact->contact_id;
                                if (!(member->final_risk_day != CB_MISSING &&
                                      t > member->final_risk_day) &&
                                    t <= member->day_exit) {
                                    sus_ptr_contact = inf_ptr_contact->pair;
                                    add_p2p_risk_history(t, member, person,
                                        sus_ptr_contact->contact_mode,
                                        sus_ptr_contact->offset,
                                        infective_prob, person->symptom);
                                }
                                inf_ptr_contact = inf_ptr_contact->next;
                            }
                        }
                    }
                }
            }
        }
    }

    /* ---- 21. Estimation ---- */
    estimation_error_type = estimation(i_inc, i_inf, 0, est, &log_likelihood,
                                       &var, &var_logit);

    /* ---- 22. Build R output ---- */
    {
        SEXP o_est, o_var, o_var_logit, o_ll, o_err;
        SEXP result, res_nm;

        PROTECT(o_est      = allocVector(REALSXP, n_par)); n_protect++;
        PROTECT(o_var      = allocVector(REALSXP, n_par * n_par)); n_protect++;
        PROTECT(o_var_logit= allocVector(REALSXP, n_par * n_par)); n_protect++;
        PROTECT(o_ll       = allocVector(REALSXP, 1)); n_protect++;
        PROTECT(o_err      = allocVector(INTSXP,  1)); n_protect++;

        for (k = 0; k < n_par; k++) REAL(o_est)[k] = (est != NULL) ? est[k] : NA_REAL;
        /* Both matrices: copy column-major for R (data[i][j] → position j*n_par+i) */
        for (i = 0; i < n_par; i++)
            for (j = 0; j < n_par; j++) {
                REAL(o_var)      [j * n_par + i] = var.data[i][j];
                REAL(o_var_logit)[j * n_par + i] = var_logit.data[i][j];
            }
        REAL(o_ll)[0]     = log_likelihood;
        INTEGER(o_err)[0] = estimation_error_type;

        PROTECT(result = allocVector(VECSXP, 5)); n_protect++;
        PROTECT(res_nm = allocVector(STRSXP, 5)); n_protect++;
        SET_VECTOR_ELT(result, 0, o_est);       SET_STRING_ELT(res_nm, 0, mkChar("est"));
        SET_VECTOR_ELT(result, 1, o_var);       SET_STRING_ELT(res_nm, 1, mkChar("var"));
        SET_VECTOR_ELT(result, 2, o_var_logit); SET_STRING_ELT(res_nm, 2, mkChar("var_logit"));
        SET_VECTOR_ELT(result, 3, o_ll);        SET_STRING_ELT(res_nm, 3, mkChar("log_likelihood"));
        SET_VECTOR_ELT(result, 4, o_err);       SET_STRING_ELT(res_nm, 4, mkChar("error_code"));
        setAttrib(result, R_NamesSymbol, res_nm);

        /* ---- 23. Cleanup ---- */
        restore_susceptibility();
        if (importance_weight != NULL) { free(importance_weight); importance_weight = NULL; }
        size_sample_states = 0;

        free_arrays();
        deflate_matrix(&var);
        deflate_matrix(&var_logit);
        deflate_matrix(&der_mat);
        if (est  != NULL) free(est);
        if (der  != NULL) free(der);
        if (der1 != NULL) free(der1);
        if (der2 != NULL) free(der2);
        if (value != NULL) free(value);
        if (p2p_covariate_buf != NULL) free(p2p_covariate_buf);
        if (pdf_incubation != NULL) { free(pdf_incubation); pdf_incubation = NULL; }
        if (sdf_incubation != NULL) { free(sdf_incubation); sdf_incubation = NULL; }

        /* free community structures */
        for (h = 0; h < n_community; h++) {
            if (community[h].contact_history != NULL) {
                for (r = 0; r < community[h].epi_duration; r++) {
                    ptr_contact = community[h].contact_history[r].c2p_contact;
                    while (ptr_contact != NULL) {
                        ptr2_contact = ptr_contact->next; free(ptr_contact);
                        ptr_contact = ptr2_contact;
                    }
                    ptr_contact = community[h].contact_history[r].p2p_contact;
                    while (ptr_contact != NULL) {
                        ptr2_contact = ptr_contact->next; free(ptr_contact);
                        ptr_contact = ptr2_contact;
                    }
                }
                free(community[h].contact_history);
            }
            if (community[h].member        != NULL) free(community[h].member);
            if (community[h].idx           != NULL) free(community[h].idx);
            if (community[h].member_impute != NULL) free(community[h].member_impute);
            if (community[h].sample_states != NULL)
                free_state_chain(community[h].sample_states, community[h].sample_states_rear);
            if (community[h].list_states   != NULL)
                free_state_chain(community[h].list_states,   community[h].list_states_rear);
            ptr_class = community[h].risk_class;
            while (ptr_class != NULL) {
                ptr_integer = ptr_class->member;
                while (ptr_integer != NULL) {
                    ptr2_integer = ptr_integer->next; free(ptr_integer);
                    ptr_integer = ptr2_integer;
                }
                if (ptr_class->risk_history != NULL) free(ptr_class->risk_history);
                ptr2_class = ptr_class->next; free(ptr_class);
                ptr_class = ptr2_class;
            }
        }
        /* free per-person arrays (must precede free(community) so that
         * community[h].epi_duration is still valid for per-person contact
         * history linked-list traversal in individualised-history mode) */
        for (i = 0; i < p_size; i++) {
            if (people[i].time_ind_covariate != NULL) free(people[i].time_ind_covariate);
            if (people[i].time_dep_covariate != NULL) free_2d_array_double(people[i].time_dep_covariate);
            if (people[i].pat_covariate      != NULL) free_2d_array_double(people[i].pat_covariate);
            if (people[i].imm_covariate      != NULL) free(people[i].imm_covariate);
            if (people[i].possible_states    != NULL) free_2d_array_int(people[i].possible_states);
            if (people[i].contact_history    != NULL) {
                h = people[i].community;
                for (r = 0; r < community[h].epi_duration; r++) {
                    ptr_contact = people[i].contact_history[r].c2p_contact;
                    while (ptr_contact != NULL) {
                        ptr2_contact = ptr_contact->next; free(ptr_contact);
                        ptr_contact = ptr2_contact;
                    }
                    ptr_contact = people[i].contact_history[r].p2p_contact;
                    while (ptr_contact != NULL) {
                        ptr2_contact = ptr_contact->next; free(ptr_contact);
                        ptr_contact = ptr2_contact;
                    }
                }
                free(people[i].contact_history);
                if (people[i].risk_history != NULL) free(people[i].risk_history);
            }
        }
        free(people); people = NULL;

        free(community); community = NULL;

        /* free cfg_pars dynamic arrays */
        if (cfg_pars.prob_incubation   != NULL) free(cfg_pars.prob_incubation);
        if (cfg_pars.prob_infectious   != NULL) free(cfg_pars.prob_infectious);
        if (cfg_pars.par_equiclass     != NULL) {
            for (k = 0; k < n_peq; k++)
                if (cfg_pars.par_equiclass[k].member != NULL)
                    free(cfg_pars.par_equiclass[k].member);
            free(cfg_pars.par_equiclass);
        }
        if (cfg_pars.ini_par_effective != NULL) free_2d_array_double(cfg_pars.ini_par_effective);
        if (cfg_pars.lower_search_bound != NULL) free(cfg_pars.lower_search_bound);
        if (cfg_pars.upper_search_bound != NULL) free(cfg_pars.upper_search_bound);
        if (cfg_pars.converge_criteria  != NULL) free(cfg_pars.converge_criteria);
        if (cfg_pars.par_fixed_id       != NULL) free(cfg_pars.par_fixed_id);
        if (cfg_pars.par_fixed_value    != NULL) free(cfg_pars.par_fixed_value);
        if (cfg_pars.c2p_covariate      != NULL) free(cfg_pars.c2p_covariate);
        if (cfg_pars.sus_p2p_covariate  != NULL) free(cfg_pars.sus_p2p_covariate);
        if (cfg_pars.inf_p2p_covariate  != NULL) free(cfg_pars.inf_p2p_covariate);
        if (cfg_pars.interaction        != NULL) free_2d_array_int(cfg_pars.interaction);
        if (cfg_pars.pat_covariate      != NULL) free(cfg_pars.pat_covariate);
        if (cfg_pars.imm_covariate      != NULL) free(cfg_pars.imm_covariate);
        if (cfg_pars.effective_lower_infectious != NULL) free(cfg_pars.effective_lower_infectious);
        if (cfg_pars.effective_upper_infectious != NULL) free(cfg_pars.effective_upper_infectious);
        if (cfg_pars.R0_multiplier      != NULL) free(cfg_pars.R0_multiplier);
        if (cfg_pars.R0_multiplier_var  != NULL) free(cfg_pars.R0_multiplier_var);

        UNPROTECT(n_protect);
        return result;
    }
}
