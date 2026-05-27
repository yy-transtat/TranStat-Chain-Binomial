#' Simulate a Single Epidemic
#'
#' Runs one chain-binomial epidemic on a pseudo-population and returns the
#' updated population list together with summary counts.  Unlike
#' \code{\link{simulate_epidemic}}, which delegates all I/O to C through
#' temporary files, this function passes the population data frames directly
#' to C and retrieves the post-simulation time-dependent covariate matrix
#' (e.g., stochastically-assigned antiviral treatment) without any file
#' intermediary.
#'
#' @param data_list Named list returned by \code{\link{gen_population}} or
#'   \code{\link{read_population}}.  Must contain \code{$pop},
#'   \code{$community}, \code{$time_ind_covariate}, and
#'   \code{$time_dep_covariate}.  Optional elements \code{$c2p_contact} and
#'   \code{$p2p_contact} are used when non-\code{NULL}; otherwise contact
#'   histories are auto-generated assuming random mixing within each community.
#' @param cfg Named list returned by \code{\link{read_config}}.  Must contain
#'   simulation parameters (\code{perform-simulation = 1},
#'   \code{simulation-only = 1}) and at least one incubation / infectious
#'   period group.
#' @param i_inc Integer (1-based).  Which incubation-period group to use.
#'   Default \code{1L}.
#' @param i_inf Integer (1-based).  Which infectious-period group to use.
#'   Default \code{1L}.
#' @param seed Integer.  RNG seed passed to the C engine.  Default
#'   \code{12345678L}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{pop}}{Updated \code{pop} data frame.  Infection, symptom,
#'     illness-onset, infection-day, index-case, and ignore fields reflect the
#'     simulated outcome.  A \code{day_infection} column is added.}
#'   \item{\code{time_ind_covariate}}{Updated time-independent covariate data
#'     frame (e.g., symptomatic-infection indicator if
#'     \code{illness-as-a-time-independent-covariate-for-infectivity} is on).}
#'   \item{\code{community}}{Updated community data frame with
#'     \code{earliest_idx_day_ill} and \code{ignore} fields reflecting
#'     case-ascertained selection.}
#'   \item{\code{time_dep_covariate}}{Updated time-dependent covariate data
#'     frame.  For studies with antiviral treatment, the treatment-status
#'     covariate is set by the simulator; this is the primary reason this
#'     function exists.}
#'   \item{\code{n_index}}{Number of index cases in the simulated epidemic.}
#'   \item{\code{n_secondary_inf}}{Number of secondary infections (non-index).}
#'   \item{\code{n_secondary_sym}}{Number of symptomatic secondary infections.}
#'   \item{\code{n_secondary_asym}}{Number of asymptomatic secondary
#'     infections.}
#'   \item{\code{n_escaped}}{Number of susceptible contacts in attacked
#'     communities who escaped infection.}
#'   \item{\code{n_preimmune}}{Number of pre-immune individuals.}
#' }
#'
#' @seealso \code{\link{gen_population}}, \code{\link{read_config}},
#'   \code{\link{simulate_epidemic}}
#'
#' @examples
#' \dontrun{
#' dat <- gen_population(
#'   n_community      = 100,
#'   community_size   = 5,
#'   day_epi_stop     = 14,
#'   case_ascertained = 1L
#' )
#' cfg <- read_config(system.file("extdata", "CaseStudy1", "config.file",
#'                                package = "ChainBinomial"))
#' out <- simulate_single(dat, cfg, seed = 42L)
#'
#' # Observed secondary attack rate
#' contacts <- subset(out$pop, idx == 0 & ignore == 0)
#' mean(contacts$infection)
#'
#' # Days of treatment (non-zero treatment covariate)
#' sum(out$time_dep_covariate$value > 0)
#' }
#'
#' @export
simulate_single <- function(data_list, cfg,
                             i_inc = 1L, i_inf = 1L,
                             seed = 12345678L) {

    stopifnot(is.list(data_list), is.list(cfg))
    stopifnot(i_inc >= 1L, i_inc <= cfg$n_inc)
    stopifnot(i_inf >= 1L, i_inf <= cfg$n_inf)

    # --- Transform sim_par_effective (raw scale → logit/log) ---
    # config.h applies: logit() for b/p/u/q parameters,
    #                   log()   for covariate coefficient parameters.
    # The boundary is: parameter index < n_b + n_p + n_u + n_q → logit.
    n_prob_pars <- cfg$n_b_mode + cfg$n_p_mode + cfg$n_u_mode + cfg$n_q_mode
    spe_tr <- numeric(cfg$n_par_equiclass)
    for (k in seq_len(cfg$n_par_equiclass)) {
        first_member_0based <- cfg$par_equiclass[[k]]$member[1L] - 1L
        raw_val <- cfg$sim_par_effective[k]
        if (first_member_0based < n_prob_pars) {
            spe_tr[k] <- log(raw_val / (1.0 - raw_val))   # logit
        } else {
            spe_tr[k] <- log(raw_val)                      # log OR
        }
    }

    # --- Transform par_fixed_value (same rule, by parameter position) ---
    fpv_tr <- numeric(0)
    if (cfg$n_par_fixed > 0L) {
        fpv_tr <- numeric(cfg$n_par_fixed)
        for (k in seq_len(cfg$n_par_fixed)) {
            pos0 <- cfg$par_fixed_id[k] - 1L
            raw  <- cfg$par_fixed_value[k]
            fpv_tr[k] <- if (pos0 < n_prob_pars) log(raw / (1.0 - raw)) else log(raw)
        }
    }

    # --- Attach transformed vectors to a local cfg copy ---
    cfg2 <- cfg
    cfg2$sim_par_effective_tr <- spe_tr
    cfg2$par_fixed_value_tr   <- if (cfg$n_par_fixed > 0L) fpv_tr else numeric(0)

    # --- Handle NA → 0 for integer fields that C expects ---
    if (is.na(cfg2$PreIllness_covariate_id)) cfg2$PreIllness_covariate_id <- 0L
    if (is.na(cfg2$Illness_covariate_id))    cfg2$Illness_covariate_id    <- 0L
    if (is.na(cfg2$RxIllness_covariate_id))  cfg2$RxIllness_covariate_id  <- 0L
    if (is.na(cfg2$RxIllness_prob))          cfg2$RxIllness_prob          <- 0.0
    if (is.na(cfg2$RxIllness_duration))      cfg2$RxIllness_duration      <- 0L
    if (is.na(cfg2$RxIllness_index_only))    cfg2$RxIllness_index_only    <- 0L

    # --- Extract and coerce contact histories from data_list ---
    # C reads the offset column with REAL(); coerce to double in case
    # read.table inferred integer for whole-number columns.
    # Shared c2p/p2p:       col 5 (community_id start stop mode offset ignore)
    # Individualised p2p:   col 6 (start stop person_i person_j mode offset ignore)
    c2p <- data_list$c2p_contact
    if (!is.null(c2p))
        c2p[[5L]] <- as.double(c2p[[5L]])

    p2p <- data_list$p2p_contact
    if (!is.null(p2p)) {
        off_col <- if (cfg$common_contact_history_within_community == 1L) 5L else 6L
        p2p[[off_col]] <- as.double(p2p[[off_col]])
    }

    # --- Call C ---
    res <- .Call("r_simulate_single",
                 data_list$pop,
                 data_list$time_ind_covariate,
                 data_list$community,
                 data_list$time_dep_covariate,
                 c2p,
                 p2p,
                 cfg2,
                 as.integer(i_inc),
                 as.integer(i_inf),
                 as.integer(seed))

    # --- Attach imputation table for asymptomatic infections ---
    res$impute <- gen_impute(res, cfg, i_inc = i_inc)

    res
}
