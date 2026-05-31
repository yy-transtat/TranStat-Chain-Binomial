#' Read a TranStat Configuration File
#'
#' Parses a TranStat \code{config.file} and returns a named list mirroring the
#' \code{CFG_PARS} C structure defined in \code{datastruct.h}.  The parsing
#' rules follow \code{config.h} exactly: each section begins with a \code{#}
#' identifier line (spaces stripped before comparison), and values are read as
#' whitespace-delimited tokens, with colons treated as whitespace.  The two
#' exceptions are \code{# input-path} and \code{# output-path}, whose values
#' are read as whole lines.
#'
#' Parameter names (the short reminder strings such as \code{b}, \code{p},
#' \code{age_sus}) that appear before numeric values in sections like
#' \code{# parameters-for-simulation} and \code{# converge-criteria} are
#' stored in companion \code{*_names} fields of the returned list.
#'
#' @param config_file Character. Path to the \code{config.file} to read.
#'
#' @return A named list with one element per \code{CFG_PARS} field that is
#'   populated by \code{config.h}.  Key fields include:
#' \describe{
#'   \item{path_in, path_out}{Input/output directory paths.}
#'   \item{n_inc, min_inc, max_inc, prob_inc}{Incubation period specification.}
#'   \item{n_inf, lower_inf, upper_inf, prob_inf}{Infectious period specification.}
#'   \item{n_b_mode, n_p_mode, n_u_mode, n_q_mode}{Numbers of transmission / group parameters.}
#'   \item{n_time_ind_covariate, n_time_dep_covariate, n_covariate}{Covariate counts.}
#'   \item{n_par_equiclass, par_equiclass}{Parameter equivalence classes (list of
#'     \code{list(size, member)}).}
#'   \item{n_par_fixed, par_fixed_id, par_fixed_value}{Fixed parameters.}
#'   \item{simulation, n_simulation, simulation_only}{Simulation switches.}
#'   \item{sim_par_names, sim_par_effective}{Simulation parameter names and values.}
#'   \item{converge_criteria_provided, converge_criteria_names, converge_criteria}{
#'     Convergence tolerances.}
#'   \item{n_ini, ini_par_provided, ini_par_names, ini_par_effective}{Initial estimates
#'     (matrix with \code{n_ini} rows and \code{n_par_equiclass} columns).}
#'   \item{search_bound_provided, search_bound_names, lower_search_bound, upper_search_bound}{
#'     Nelder-Mead search bounds.}
#'   \item{EM}{Whether the EM algorithm is used (1) or not (0).}
#' }
#'
#' @examples
#' cfg_file <- system.file("extdata", "CaseStudy1", "config.file", package = "ChainBinomial")
#' cfg <- read_config(cfg_file)
#' cfg$n_b_mode        # number of c2p transmission probabilities
#' cfg$sim_par_names   # parameter names used in simulation
#' cfg$sim_par_effective
#'
#' @export
read_config <- function(config_file) {
    lines <- readLines(config_file, warn = FALSE)

    # -----------------------------------------------------------------------
    # 1.  Split file into named sections
    #     Each line that starts with '#' opens a new section.
    #     The identifier is the text after '#' with all spaces removed
    #     (matching RemoveSpaces() in config.h).
    # -----------------------------------------------------------------------
    is_hdr  <- grepl("^\\s*#", lines)
    starts  <- which(is_hdr)
    ends    <- c(starts[-1] - 1L, length(lines))

    sec_ids <- vapply(starts, function(i)
        gsub("\\s+", "", sub("^\\s*#", "", lines[i])), character(1))

    sec_bodies <- lapply(seq_along(starts), function(k) {
        i1 <- starts[k] + 1L
        i2 <- ends[k]
        if (i1 <= i2) lines[i1:i2] else character(0)
    })
    names(sec_bodies) <- sec_ids

    # -----------------------------------------------------------------------
    # 2.  Helper factories
    # -----------------------------------------------------------------------

    # Return the body lines of a section (empty vector if absent).
    body <- function(id) sec_bodies[[id]] %||% character(0)
    `%||%` <- function(x, y) if (!is.null(x)) x else y

    # First non-empty line (used for whole-line path reads).
    first_line <- function(id) {
        for (l in body(id)) if (nchar(trimws(l)) > 0L) return(trimws(l))
        ""
    }

    # Build a token stream from a section.
    # Colons are treated as whitespace (mirrors fscanf "%d:" behaviour).
    tok_stream <- function(id) {
        raw  <- paste(body(id), collapse = " ")
        raw  <- gsub(":", " ", raw)
        toks <- strsplit(trimws(raw), "\\s+")[[1]]
        toks <- toks[nchar(toks) > 0L]
        pos  <- 1L
        list(
            gi  = function() { v <- suppressWarnings(as.integer(toks[pos])); pos <<- pos + 1L; v },
            gd  = function() { v <- suppressWarnings(as.double(toks[pos]));  pos <<- pos + 1L; v },
            gs  = function() { v <- toks[pos];                               pos <<- pos + 1L; v },
            ok  = function() pos <= length(toks)
        )
    }

    # Read a single integer / double from a section (returns 0 if absent).
    gi1 <- function(id) { s <- tok_stream(id); if (s$ok()) s$gi() else 0L }
    gd1 <- function(id) { s <- tok_stream(id); if (s$ok()) s$gd() else 0.0 }

    # Parse "n: id1 id2 ..." covariate-index sections.
    parse_cov_ids <- function(id) {
        s <- tok_stream(id)
        n <- if (s$ok()) s$gi() else 0L
        ids <- if (!is.na(n) && n > 0L)
            vapply(seq_len(n), function(.) s$gi(), integer(1))
        else integer(0)
        list(n = n %||% 0L, ids = ids)
    }

    cfg <- list()

    # -----------------------------------------------------------------------
    # 3.  Paths  (whole-line reads)
    # -----------------------------------------------------------------------
    cfg$path_in  <- first_line("input-path")
    cfg$path_out <- first_line("output-path")

    # -----------------------------------------------------------------------
    # 4.  Incubation period
    #     Format: n_inc
    #             min  max  (prob_0  prob_1 ... prob_{len-1})  <-- repeat n_inc times
    # -----------------------------------------------------------------------
    s            <- tok_stream("min-max-days-and-probs-of-incubation-period")
    cfg$n_inc    <- s$gi()
    cfg$min_inc  <- integer(cfg$n_inc)
    cfg$max_inc  <- integer(cfg$n_inc)
    cfg$prob_inc <- vector("list", cfg$n_inc)
    for (i in seq_len(cfg$n_inc)) {
        cfg$min_inc[i] <- s$gi()
        cfg$max_inc[i] <- s$gi()
        len <- cfg$max_inc[i] - cfg$min_inc[i] + 1L
        cfg$prob_inc[[i]] <- vapply(seq_len(len), function(.) s$gd(), double(1))
    }

    # -----------------------------------------------------------------------
    # 5.  Infectious period
    #     Format: n_inf
    #             lower  upper  (prob_0 ... prob_{len-1})  <-- repeat n_inf times
    # -----------------------------------------------------------------------
    s             <- tok_stream("primary-lower-upper-bounds-and-probs-of-infectious-days-relative-to-symptom-onset-day")
    cfg$n_inf     <- s$gi()
    cfg$lower_inf <- integer(cfg$n_inf)
    cfg$upper_inf <- integer(cfg$n_inf)
    cfg$prob_inf  <- vector("list", cfg$n_inf)
    for (i in seq_len(cfg$n_inf)) {
        cfg$lower_inf[i] <- s$gi()
        cfg$upper_inf[i] <- s$gi()
        len <- cfg$upper_inf[i] - cfg$lower_inf[i] + 1L
        cfg$prob_inf[[i]] <- vapply(seq_len(len), function(.) s$gd(), double(1))
    }

    # -----------------------------------------------------------------------
    # 6.  Basic mode / group counts
    # -----------------------------------------------------------------------
    cfg$n_b_mode             <- gi1("number-of-c2p-transmission-probabilities")
    cfg$n_p_mode             <- gi1("number-of-p2p-transmission-probabilities")
    cfg$n_u_mode             <- gi1("number-of-pathogenicity-groups")
    cfg$n_q_mode             <- gi1("number-of-preseason-immunity-groups")
    cfg$n_time_ind_covariate <- gi1("number-of-time-independent-covariates")
    cfg$n_time_dep_covariate <- gi1("number-of-time-dependent-covariates")
    cfg$n_covariate          <- cfg$n_time_ind_covariate + cfg$n_time_dep_covariate

    # -----------------------------------------------------------------------
    # 7.  Special covariate flags
    # -----------------------------------------------------------------------
    # Pre/post illness onset as time-dependent covariate for infectivity
    cfg$PreIllness_as_covariate <- 0L
    cfg$PreIllness_covariate_id <- NA_integer_
    s <- tok_stream("pre-post-illness-onset-as-a-time-dependent-covariate-for-infectivity")
    if (s$ok()) {
        cfg$PreIllness_as_covariate <- s$gi()
        if (isTRUE(cfg$PreIllness_as_covariate == 1L) && s$ok())
            cfg$PreIllness_covariate_id <- s$gi()
    }

    # Illness as time-independent covariate for infectivity
    cfg$Illness_as_covariate <- 0L
    cfg$Illness_covariate_id <- NA_integer_
    s <- tok_stream("illness-as-a-time-independent-covariate-for-infectivity")
    if (s$ok()) {
        cfg$Illness_as_covariate <- s$gi()
        if (isTRUE(cfg$Illness_as_covariate == 1L) && s$ok())
            cfg$Illness_covariate_id <- s$gi()
    }

    # Antiviral treatment as time-dependent covariate
    # Format: as_covariate: [covariate_id  prob  duration  index_only]
    cfg$RxIllness_as_covariate <- 0L
    cfg$RxIllness_covariate_id <- NA_integer_
    cfg$RxIllness_prob         <- NA_real_
    cfg$RxIllness_duration     <- NA_integer_
    cfg$RxIllness_index_only   <- NA_integer_
    s <- tok_stream("antiviral-treatment-as-a-time-dependent-covariate-for-infectivity")
    if (s$ok()) {
        cfg$RxIllness_as_covariate <- s$gi()
        if (isTRUE(cfg$RxIllness_as_covariate == 1L)) {
            cfg$RxIllness_covariate_id <- s$gi()
            cfg$RxIllness_prob         <- s$gd()
            cfg$RxIllness_duration     <- s$gi()
            cfg$RxIllness_index_only   <- s$gi()
        }
    }

    # -----------------------------------------------------------------------
    # 8.  Covariate index arrays
    #     Each section format: n: id_1 id_2 ... id_n
    # -----------------------------------------------------------------------
    tmp <- parse_cov_ids("covariates-affecting-susceptibility-for-c2p-transmission")
    cfg$n_c2p_covariate  <- tmp$n;  cfg$c2p_covariate      <- tmp$ids

    tmp <- parse_cov_ids("covariates-affecting-susceptibility-for-p2p-transmission")
    cfg$n_sus_p2p_covariate <- tmp$n; cfg$sus_p2p_covariate <- tmp$ids

    tmp <- parse_cov_ids("covariates-affecting-infectiousness-for-p2p-transmission")
    cfg$n_inf_p2p_covariate <- tmp$n; cfg$inf_p2p_covariate <- tmp$ids

    # Interactions: n: (sus_idx inf_idx) pairs
    s <- tok_stream("interactions-for-p2p-transmission")
    cfg$n_int_p2p_covariate <- if (s$ok()) s$gi() else 0L
    cfg$interaction <- if (isTRUE(cfg$n_int_p2p_covariate > 0L))
        lapply(seq_len(cfg$n_int_p2p_covariate), function(.) c(s$gi(), s$gi()))
    else NULL

    tmp <- parse_cov_ids("covariates-affecting-pathogenicity")
    cfg$n_pat_covariate <- tmp$n; cfg$pat_covariate <- tmp$ids

    tmp <- parse_cov_ids("covariates-affecting-preseason-immunity")
    cfg$n_imm_covariate <- tmp$n; cfg$imm_covariate <- tmp$ids

    # Derived counts (same formulas as config.h)
    cfg$n_p2p_covariate <- cfg$n_sus_p2p_covariate + cfg$n_inf_p2p_covariate + cfg$n_int_p2p_covariate
    cfg$n_par <- cfg$n_b_mode + cfg$n_p_mode + cfg$n_u_mode + cfg$n_q_mode +
                 cfg$n_c2p_covariate + cfg$n_p2p_covariate +
                 cfg$n_pat_covariate + cfg$n_imm_covariate

    # -----------------------------------------------------------------------
    # 9.  Parameter equivalence classes
    #     Format: n_equiclass:
    #             size_0: member_0_0 [member_0_1 ...]
    #             ...
    # -----------------------------------------------------------------------
    s <- tok_stream("equal-parameters")
    cfg$n_par_equiclass <- if (s$ok()) s$gi() else 0L
    cfg$par_equiclass <- if (isTRUE(cfg$n_par_equiclass > 0L)) {
        lapply(seq_len(cfg$n_par_equiclass), function(.) {
            sz <- s$gi()
            list(size   = sz,
                 member = vapply(seq_len(sz), function(.) s$gi(), integer(1)))
        })
    } else NULL

    # -----------------------------------------------------------------------
    # 10.  Fixed parameters
    #      Format: n: (id: value) pairs
    # -----------------------------------------------------------------------
    s <- tok_stream("fixed-parameters")
    cfg$n_par_fixed     <- if (s$ok()) s$gi() else 0L
    cfg$par_fixed_id    <- integer(0)
    cfg$par_fixed_value <- double(0)
    if (isTRUE(cfg$n_par_fixed > 0L)) {
        cfg$par_fixed_id    <- integer(cfg$n_par_fixed)
        cfg$par_fixed_value <- double(cfg$n_par_fixed)
        for (i in seq_len(cfg$n_par_fixed)) {
            cfg$par_fixed_id[i]    <- s$gi()
            cfg$par_fixed_value[i] <- s$gd()
        }
    }

    # -----------------------------------------------------------------------
    # 11.  Simulation parameters
    # -----------------------------------------------------------------------
    cfg$simulation                                <- gi1("perform-simulation")
    cfg$n_simulation                              <- gi1("number-of-simulations")
    cfg$simulation_only                           <- gi1("simulation-only")
    cfg$output_simulation_data                    <- gi1("output-simulation-data")
    cfg$idx_initiated_followup_in_simulation      <- gi1("index-case-initiate-followup-in-simulation")
    cfg$followup_duration_after_idx_in_simulation <- gi1("followup-duration-after-index-case-in-simulation")
    cfg$prop_mix_imm_esc                          <- gd1("proportion-with-ambiguity-about-preimmunity-and-escape-status")
    cfg$asym_effect_sim                           <- gd1("relative-infectivity-of-asymptomatic-case-for-simulation")

    # Simulation parameter values
    # Format: name_0  value_0 / name_1  value_1 / ...  (one pair per equivalence class)
    cfg$sim_par_names     <- character(0)
    cfg$sim_par_effective <- double(0)
    s <- tok_stream("parameters-for-simulation")
    if (s$ok() && isTRUE(cfg$n_par_equiclass > 0L)) {
        cfg$sim_par_names     <- character(cfg$n_par_equiclass)
        cfg$sim_par_effective <- double(cfg$n_par_equiclass)
        for (i in seq_len(cfg$n_par_equiclass)) {
            cfg$sim_par_names[i]     <- s$gs()
            cfg$sim_par_effective[i] <- s$gd()
        }
    }

    # -----------------------------------------------------------------------
    # 12.  Estimation parameters
    # -----------------------------------------------------------------------
    cfg$optimization_choice <- gi1("optimization-choice")

    # Convergence criteria
    # Format: provided:
    #         name_0  tol_0 / ...  (one pair per equivalence class)
    s <- tok_stream("converge-criteria")
    cfg$converge_criteria_provided <- if (s$ok()) s$gi() else 0L
    cfg$converge_criteria_names    <- character(cfg$n_par_equiclass)
    cfg$converge_criteria          <- double(cfg$n_par_equiclass)
    if (isTRUE(cfg$converge_criteria_provided == 1L) && isTRUE(cfg$n_par_equiclass > 0L)) {
        for (i in seq_len(cfg$n_par_equiclass)) {
            cfg$converge_criteria_names[i] <- s$gs()
            cfg$converge_criteria[i]       <- s$gd()
        }
    }

    # Initial estimates
    # Format: n_ini:ini_par_provided
    #         name_0  val_0 / ...  (repeated n_ini times, one row per set)
    s <- tok_stream("initial-estimates")
    cfg$n_ini            <- if (s$ok()) s$gi() else 1L
    cfg$ini_par_provided <- if (s$ok()) s$gi() else 0L
    cfg$ini_par_names    <- character(cfg$n_par_equiclass)
    cfg$ini_par_effective <- matrix(0.0, nrow = cfg$n_ini, ncol = cfg$n_par_equiclass)
    if (isTRUE(cfg$ini_par_provided == 1L) && isTRUE(cfg$n_par_equiclass > 0L)) {
        for (i in seq_len(cfg$n_ini)) {
            for (j in seq_len(cfg$n_par_equiclass)) {
                nm <- s$gs()
                if (i == 1L) cfg$ini_par_names[j] <- nm
                cfg$ini_par_effective[i, j] <- s$gd()
            }
        }
    }

    # Search bounds (for Nelder-Mead)
    # Format: provided:
    #         name_0  lower_0  upper_0 / ...
    s <- tok_stream("search-bounds")
    cfg$search_bound_provided <- if (s$ok()) s$gi() else 0L
    cfg$search_bound_names    <- character(cfg$n_par_equiclass)
    cfg$lower_search_bound    <- double(cfg$n_par_equiclass)
    cfg$upper_search_bound    <- double(cfg$n_par_equiclass)
    if (isTRUE(cfg$search_bound_provided == 1L) && isTRUE(cfg$n_par_equiclass > 0L)) {
        for (i in seq_len(cfg$n_par_equiclass)) {
            cfg$search_bound_names[i] <- s$gs()
            cfg$lower_search_bound[i] <- s$gd()
            cfg$upper_search_bound[i] <- s$gd()
        }
    }

    # EM algorithm
    cfg$EM                           <- gi1("perform-EM-algorithm")
    cfg$min_size_MCEM                <- gi1("min-number-of-possible-status-to-use-mcem")
    cfg$community_specific_weighting <- gi1("use-community-specific-weighting")
    cfg$n_base_sampling              <- gi1("number-of-base-mcmc-samples")
    cfg$n_burnin_sampling            <- gi1("number-of-burnin-mcmc-samples")
    cfg$n_burnin_iter                <- gi1("number-of-burnin-mcmc-iterations")
    cfg$n_sampling_for_mce           <- gi1("number-of-samplings-for-mc-error")
    cfg$use_bootstrap_for_mce        <- gi1("use-bootstrap-for-mc-error")
    cfg$skip_Evar_for_mce            <- gi1("do-not-calculate-average-variance-for-mc-error")
    cfg$check_missingness            <- gi1("check-missingness")
    cfg$check_mixing                 <- gi1("check-mixing")
    cfg$check_runtime                <- gi1("check-runtime")

    cfg$asym_effect_est                         <- gd1("relative-infectivity-of-asymptomatic-case-for-estimation")
    cfg$common_contact_history_within_community <- gi1("members-share-common-contact-history-within-communities")
    cfg$generate_c2p_contact                    <- gi1("automatically-generate-c2p-contact-file")
    cfg$generate_p2p_contact                    <- gi1("automatically-generate-p2p-contact-file")
    cfg$c2p_offset                              <- gi1("use-c2p-offset")
    cfg$p2p_offset                              <- gi1("use-p2p-offset")
    cfg$adjust_for_left_truncation              <- gi1("adjust-for-selection-bias")
    cfg$adjust_for_right_censoring              <- gi1("adjust-for-right-censoring")
    cfg$preset_index                            <- gi1("prefix-index-cases")

    # Use index case onset days to improve estimation of b
    # Format: use_flag: [n_c2p_group]
    s <- tok_stream("use-onset-days-of-index-cases-to-improve-estimation-of-b")
    cfg$use_index_cases_to_improve_b <- if (s$ok()) s$gi() else 0L
    cfg$n_c2p_group <- 0L
    cfg$c2p_group   <- NULL
    if (isTRUE(cfg$use_index_cases_to_improve_b == 1L) && s$ok())
        cfg$n_c2p_group <- s$gi()

    cfg$CPI_duration <- gi1("epidemic-duration-for-calculating-CPI")

    # Effective lower/upper infectious bounds (for SAR/R0 calculation)
    # Format: provided: [lower_0 upper_0 ...]  (one pair per p_mode)
    s <- tok_stream("effective-lower-upper-bounds-of-infectious-days-relative-to-symptom-onset-day")
    cfg$effective_bounds_provided  <- if (s$ok()) s$gi() else 0L
    cfg$effective_lower_infectious <- NULL
    cfg$effective_upper_infectious <- NULL
    if (isTRUE(cfg$effective_bounds_provided == 1L) && isTRUE(cfg$n_p_mode > 0L)) {
        cfg$effective_lower_infectious <- integer(cfg$n_p_mode)
        cfg$effective_upper_infectious <- integer(cfg$n_p_mode)
        for (i in seq_len(cfg$n_p_mode)) {
            cfg$effective_lower_infectious[i] <- s$gi()
            cfg$effective_upper_infectious[i] <- s$gi()
        }
    }

    # -----------------------------------------------------------------------
    # 13.  SAR covariate sets
    #      Format: n_sets:
    #        [if n_time_ind_covariate > 0, for each set:]
    #          sus-ind  v0 v1 ... v{n_tic-1}
    #          inf-ind  v0 v1 ... v{n_tic-1}
    #        [if n_time_dep_covariate > 0:]
    #          lower  upper
    #          [for each set:]
    #            sus-dep  m  {m × rows of: start stop v0..v{n_tdc-1}}
    #            inf-dep  m  {m × rows of: start stop v0..v{n_tdc-1}}
    # -----------------------------------------------------------------------

    s <- tok_stream("covariates-for-calculating-SAR-provided")

    cfg$SAR_n_covariate_sets       <- if (s$ok()) s$gi() else 0L
    cfg$SAR_sus_time_ind_covariate <- NULL
    cfg$SAR_inf_time_ind_covariate <- NULL
    cfg$SAR_time_dep_lower         <- 0L
    cfg$SAR_time_dep_upper         <- 0L
    cfg$SAR_sus_time_dep_covariate <- NULL
    cfg$SAR_inf_time_dep_covariate <- NULL

    if (isTRUE(cfg$SAR_n_covariate_sets > 0L) && isTRUE(cfg$n_p_mode > 0L)) {
        n_sets <- cfg$SAR_n_covariate_sets
        n_tic  <- cfg$n_time_ind_covariate
        n_tdc  <- cfg$n_time_dep_covariate

        # ---- time-independent covariates (one sus-ind / inf-ind block per set) ----
        if (n_tic > 0L) {
            sti <- matrix(0.0, nrow = n_sets, ncol = n_tic)
            iti <- matrix(0.0, nrow = n_sets, ncol = n_tic)
            for (nn in seq_len(n_sets)) {
                s$gs()   # consume "sus-ind" label
                for (i in seq_len(n_tic)) sti[nn, i] <- s$gd()
                s$gs()   # consume "inf-ind" label
                for (i in seq_len(n_tic)) iti[nn, i] <- s$gd()
            }
            cfg$SAR_sus_time_ind_covariate <- sti
            cfg$SAR_inf_time_ind_covariate <- iti
        }

        # ---- time-dependent covariates ----
        if (n_tdc > 0L) {
            cfg$SAR_time_dep_lower <- s$gi()
            cfg$SAR_time_dep_upper <- s$gi()
            len <- cfg$SAR_time_dep_upper - cfg$SAR_time_dep_lower + 1L

            std <- vector("list", n_sets)   # [[set]][time_row, cov_col]
            itd <- vector("list", n_sets)
            for (nn in seq_len(n_sets)) {
                std[[nn]] <- matrix(0.0, nrow = len, ncol = n_tdc)
                itd[[nn]] <- matrix(0.0, nrow = len, ncol = n_tdc)
            }

            for (nn in seq_len(n_sets)) {
                # susceptible side
                s$gs()              # "sus-dep"
                m_sus <- s$gi()
                for (i in seq_len(m_sus)) {
                    st   <- s$gi(); sp <- s$gi()
                    vals <- vapply(seq_len(n_tdc), function(.) s$gd(), double(1L))
                    if (sp >= st) {
                        for (t in seq.int(st, sp)) {
                            r <- t - cfg$SAR_time_dep_lower + 1L
                            std[[nn]][r, ] <- vals
                        }
                    }
                }
                # infectious side
                s$gs()              # "inf-dep"
                m_inf <- s$gi()
                for (i in seq_len(m_inf)) {
                    st   <- s$gi(); sp <- s$gi()
                    vals <- vapply(seq_len(n_tdc), function(.) s$gd(), double(1L))
                    if (sp >= st) {
                        for (t in seq.int(st, sp)) {
                            r <- t - cfg$SAR_time_dep_lower + 1L
                            itd[[nn]][r, ] <- vals
                        }
                    }
                }
            }

            cfg$SAR_sus_time_dep_covariate <- std
            cfg$SAR_inf_time_dep_covariate <- itd
        }
    }

    # -----------------------------------------------------------------------
    # 14.  R0 multiplier
    #      Format: provided: [mult_0 var_0 ...]  (one pair per p_mode)
    # -----------------------------------------------------------------------
    s <- tok_stream("multiplier-for-calculating-R0")
    cfg$R0_multiplier_provided <- if (s$ok()) s$gi() else 0L
    cfg$R0_multiplier          <- NULL
    cfg$R0_multiplier_var      <- NULL
    if (isTRUE(cfg$R0_multiplier_provided > 0L) && isTRUE(cfg$n_p_mode > 0L)) {
        cfg$R0_multiplier     <- double(cfg$n_p_mode)
        cfg$R0_multiplier_var <- double(cfg$n_p_mode)
        for (i in seq_len(cfg$n_p_mode)) {
            cfg$R0_multiplier[i]     <- s$gd()
            cfg$R0_multiplier_var[i] <- s$gd()
        }
    }

    # Serial division for time-varying R0
    # Format: divide_flag: [start stop]  or  [start stop window]
    s <- tok_stream("serial-division-of-epidemic-for-calculating-time-varying-R0")
    cfg$R0_divide_by_time <- if (s$ok()) s$gi() else 0L
    cfg$R0_divide_start   <- 0L
    cfg$R0_divide_stop    <- 0L
    cfg$R0_window_size    <- 0L
    if (isTRUE(cfg$R0_divide_by_time == 1L) && isTRUE(cfg$n_p_mode > 0L)) {
        cfg$R0_divide_start <- s$gi()
        cfg$R0_divide_stop  <- s$gi()
        if (cfg$n_p_mode != 2L) cfg$R0_window_size <- s$gi()
    }

    # -----------------------------------------------------------------------
    # 15.  Output / diagnostic switches
    # -----------------------------------------------------------------------
    cfg$goodness_of_fit     <- gi1("goodness-of-fit")
    cfg$stat_test           <- gi1("perform-statistical-test")
    cfg$skip_variance       <- gi1("do-not-estimate-variance")
    cfg$print_covariance    <- gi1("print-covariance-matrix")
    cfg$skip_output         <- gi1("do-not-output-estimates")
    cfg$simplify_output     <- gi1("simplify-output")
    cfg$simplify_output_SAR <- gi1("simplify-output-SAR")
    cfg$simplify_output_R0  <- gi1("simplify-output-R0")
    cfg$silent_run          <- gi1("run-transtat-silently")
    cfg$write_error_log     <- gi1("write-error-log")

    # -----------------------------------------------------------------------
    # 16.  Covariate labels and parameter labels  (R-only; not read by C)
    #
    #      Delegates to update_var_par_labels() with data_list = NULL so that
    #      placeholder labels ("x1", "x2", ..., "c2p_x1", "p2p_i_x2", ...)
    #      are set immediately.  Call update_var_par_labels(cfg, data_list)
    #      after read_population() to replace placeholders with real names.
    # -----------------------------------------------------------------------
    cfg <- update_var_par_labels(cfg)

    cfg
}
