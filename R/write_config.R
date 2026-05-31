#' Write a Configuration Object to a TranStat Config File
#'
#' Serialises a configuration list (as returned by \code{\link{read_config}} or
#' \code{\link{set_config}}) to a text file in the format expected by the
#' TranStat C engine.  The written file can be read back verbatim with
#' \code{\link{read_config}}.
#'
#' @param cfg Named list returned by \code{\link{read_config}} or
#'   \code{\link{set_config}}.
#' @param dir Character. Directory in which to write the file.  Created
#'   recursively if it does not exist.
#' @param file_name Character. Name of the output file (e.g.\
#'   \code{"config.file"}).  \strong{Must} be supplied; \code{NULL} throws an
#'   error to prevent accidentally overwriting an existing configuration file.
#'
#' @return The full path to the written file, invisibly.
#'
#' @seealso \code{\link{read_config}}, \code{\link{set_config}}
#'
#' @examples
#' \dontrun{
#' cfg  <- read_config(system.file("extdata", "CaseStudy1", "config.file",
#'                                 package = "ChainBinomial"))
#' cfg2 <- set_config(cfg, n_base_sampling = 1000L)
#' write_config(cfg2, dir = tempdir(), file_name = "config_new.file")
#' }
#'
#' @export
write_config <- function(cfg, dir, file_name) {

    if (is.null(file_name))
        stop("file_name must be provided; no default is set to avoid ",
             "accidentally overwriting an existing configuration file.")
    stopifnot(is.list(cfg), is.character(dir), is.character(file_name))

    `%||%` <- function(x, y) if (!is.null(x)) x else y

    if (!dir.exists(dir))
        dir.create(dir, recursive = TRUE)

    out_path <- file.path(dir, file_name)
    lns      <- character(0)

    # ---- internal helpers --------------------------------------------------

    # Append section: "# header\n<body_lines>\n"
    sec <- function(header, body_lines) {
        lns <<- c(lns, paste0("# ", header))
        if (length(body_lines) > 0L)
            lns <<- c(lns, body_lines)
        lns <<- c(lns, "")
    }

    fi  <- function(x) as.character(as.integer(x))
    fd  <- function(x) sprintf("%.10g", as.double(x))
    fvi <- function(x) paste(fi(x), collapse = " ")
    fvd <- function(x) paste(fd(x), collapse = " ")

    # "n: id1 id2 ..." covariate-index section (or "0:")
    sec_cov <- function(header, n, ids) {
        n <- as.integer(n %||% 0L)
        sec(header, if (n > 0L) paste0(fi(n), ":", fvi(ids)) else "0:")
    }

    # ---- 1.  Paths ---------------------------------------------------------
    sec("input-path",  cfg$path_in  %||% "")
    sec("output-path", cfg$path_out %||% "")

    # ---- 2.  Incubation period ---------------------------------------------
    {
        body <- fi(cfg$n_inc)
        for (i in seq_len(cfg$n_inc)) {
            body <- c(body,
                      paste(fi(cfg$min_inc[i]), fi(cfg$max_inc[i])),
                      fvd(cfg$prob_inc[[i]]),
                      "")
        }
        sec("min-max-days-and-probs-of-incubation-period", body)
    }

    # ---- 3.  Infectious period ---------------------------------------------
    {
        body <- fi(cfg$n_inf)
        for (i in seq_len(cfg$n_inf)) {
            body <- c(body,
                      paste(fi(cfg$lower_inf[i]), fi(cfg$upper_inf[i])),
                      fvd(cfg$prob_inf[[i]]),
                      "")
        }
        sec("primary-lower-upper-bounds-and-probs-of-infectious-days-relative-to-symptom-onset-day",
            body)
    }

    # ---- 4.  Mode / group counts -------------------------------------------
    sec("number-of-c2p-transmission-probabilities", fi(cfg$n_b_mode             %||% 0L))
    sec("number-of-p2p-transmission-probabilities", fi(cfg$n_p_mode             %||% 0L))
    sec("number-of-pathogenicity-groups",           fi(cfg$n_u_mode             %||% 0L))
    sec("number-of-preseason-immunity-groups",      fi(cfg$n_q_mode             %||% 0L))
    sec("number-of-time-independent-covariates",    fi(cfg$n_time_ind_covariate %||% 0L))
    sec("number-of-time-dependent-covariates",      fi(cfg$n_time_dep_covariate %||% 0L))

    # ---- 5.  Covariate index arrays ----------------------------------------
    sec_cov("covariates-affecting-susceptibility-for-c2p-transmission",
            cfg$n_c2p_covariate, cfg$c2p_covariate)
    sec_cov("covariates-affecting-susceptibility-for-p2p-transmission",
            cfg$n_sus_p2p_covariate, cfg$sus_p2p_covariate)
    sec_cov("covariates-affecting-infectiousness-for-p2p-transmission",
            cfg$n_inf_p2p_covariate, cfg$inf_p2p_covariate)

    {
        n <- as.integer(cfg$n_int_p2p_covariate %||% 0L)
        if (n > 0L) {
            pairs_str <- paste(
                vapply(seq_len(n),
                       function(j) paste(fi(cfg$interaction[[j]]), collapse = " "),
                       character(1L)),
                collapse = " ")
            sec("interactions-for-p2p-transmission",
                paste0(fi(n), ":", pairs_str))
        } else {
            sec("interactions-for-p2p-transmission", "0:")
        }
    }

    sec_cov("covariates-affecting-pathogenicity",
            cfg$n_pat_covariate, cfg$pat_covariate)
    sec_cov("covariates-affecting-preseason-immunity",
            cfg$n_imm_covariate, cfg$imm_covariate)

    # ---- 6.  Special covariate flags ---------------------------------------

    # Illness as TIC for infectivity
    {
        ill <- as.integer(cfg$Illness_as_covariate %||% 0L)
        if (isTRUE(ill == 1L) && !is.na(cfg$Illness_covariate_id %||% NA))
            sec("illness-as-a-time-independent-covariate-for-infectivity",
                paste0(fi(ill), ":", fi(cfg$Illness_covariate_id)))
        else
            sec("illness-as-a-time-independent-covariate-for-infectivity",
                paste0(fi(ill), ":"))
    }

    # Pre/post illness onset as TDC for infectivity
    {
        pre <- as.integer(cfg$PreIllness_as_covariate %||% 0L)
        if (isTRUE(pre == 1L) && !is.na(cfg$PreIllness_covariate_id %||% NA))
            sec("pre-post-illness-onset-as-a-time-dependent-covariate-for-infectivity",
                paste0(fi(pre), ":", fi(cfg$PreIllness_covariate_id)))
        else
            sec("pre-post-illness-onset-as-a-time-dependent-covariate-for-infectivity",
                paste0(fi(pre), ":"))
    }

    # Antiviral treatment as TDC for infectivity
    {
        rx <- as.integer(cfg$RxIllness_as_covariate %||% 0L)
        if (isTRUE(rx == 1L))
            sec("antiviral-treatment-as-a-time-dependent-covariate-for-infectivity",
                paste(fi(rx), fi(cfg$RxIllness_covariate_id),
                      fd(cfg$RxIllness_prob),
                      fi(cfg$RxIllness_duration),
                      fi(cfg$RxIllness_index_only)))
        else
            sec("antiviral-treatment-as-a-time-dependent-covariate-for-infectivity",
                paste0(fi(rx), ":"))
    }

    # ---- 7.  Parameter equivalence classes ---------------------------------
    {
        n_ec <- as.integer(cfg$n_par_equiclass %||% 0L)
        if (n_ec > 0L) {
            body <- paste0(fi(n_ec), ":")
            for (k in seq_len(n_ec)) {
                ec   <- cfg$par_equiclass[[k]]
                body <- c(body, paste0(fi(ec$size), ":", fvi(ec$member)))
            }
            sec("equal-parameters", body)
        } else {
            sec("equal-parameters", "0:")
        }
    }

    # ---- 8.  Fixed parameters ----------------------------------------------
    {
        n_fp <- as.integer(cfg$n_par_fixed %||% 0L)
        if (n_fp > 0L) {
            body <- paste0(fi(n_fp), ":")
            for (k in seq_len(n_fp))
                body <- c(body, paste(fi(cfg$par_fixed_id[k]),
                                      fd(cfg$par_fixed_value[k])))
            sec("fixed-parameters", body)
        } else {
            sec("fixed-parameters", "0:")
        }
    }

    # ---- 9.  Simulation parameters -----------------------------------------
    sec("perform-simulation",     fi(cfg$simulation             %||% 0L))
    sec("output-simulation-data", fi(cfg$output_simulation_data %||% 0L))
    sec("number-of-simulations",  fi(cfg$n_simulation           %||% 0L))
    sec("simulation-only",        fi(cfg$simulation_only        %||% 0L))

    {
        n_ec <- as.integer(cfg$n_par_equiclass %||% 0L)
        if (n_ec > 0L && length(cfg$sim_par_names) == n_ec)
            sec("parameters-for-simulation",
                vapply(seq_len(n_ec),
                       function(j) paste(cfg$sim_par_names[j],
                                         fd(cfg$sim_par_effective[j])),
                       character(1L)))
        else
            sec("parameters-for-simulation", character(0L))
    }

    sec("index-case-initiate-followup-in-simulation",
        fi(cfg$idx_initiated_followup_in_simulation      %||% 0L))
    sec("followup-duration-after-index-case-in-simulation",
        fi(cfg$followup_duration_after_idx_in_simulation %||% 0L))
    sec("proportion-with-ambiguity-about-preimmunity-and-escape-status",
        fd(cfg$prop_mix_imm_esc %||% 0))

    # ---- 10.  Estimation parameters ----------------------------------------
    sec("optimization-choice", fi(cfg$optimization_choice %||% 0L))

    {
        prov <- as.integer(cfg$converge_criteria_provided %||% 0L)
        n_ec <- as.integer(cfg$n_par_equiclass %||% 0L)
        body <- paste0(fi(prov), ":")
        if (isTRUE(prov == 1L) && n_ec > 0L)
            body <- c(body, vapply(seq_len(n_ec), function(j)
                paste(cfg$converge_criteria_names[j],
                      fd(cfg$converge_criteria[j])),
                character(1L)))
        sec("converge-criteria", body)
    }

    {
        n_ini <- as.integer(cfg$n_ini           %||% 1L)
        prov  <- as.integer(cfg$ini_par_provided %||% 0L)
        n_ec  <- as.integer(cfg$n_par_equiclass  %||% 0L)
        body  <- paste0(fi(n_ini), ":", fi(prov))
        if (isTRUE(prov == 1L) && n_ec > 0L) {
            for (i in seq_len(n_ini)) {
                body <- c(body,
                          vapply(seq_len(n_ec),
                                 function(j) paste(cfg$ini_par_names[j],
                                                   fd(cfg$ini_par_effective[i, j])),
                                 character(1L)))
                if (i < n_ini) body <- c(body, "")
            }
        }
        sec("initial-estimates", body)
    }

    {
        prov <- as.integer(cfg$search_bound_provided %||% 0L)
        n_ec <- as.integer(cfg$n_par_equiclass %||% 0L)
        body <- paste0(fi(prov), ":")
        if (isTRUE(prov == 1L) && n_ec > 0L)
            body <- c(body, vapply(seq_len(n_ec), function(j)
                paste(cfg$search_bound_names[j],
                      fd(cfg$lower_search_bound[j]),
                      fd(cfg$upper_search_bound[j])),
                character(1L)))
        sec("search-bounds", body)
    }

    # ---- 11.  EM / MCMC parameters -----------------------------------------
    sec("perform-EM-algorithm",
        fi(cfg$EM %||% 1L))
    sec("min-number-of-possible-status-to-use-mcem",
        fi(cfg$min_size_MCEM %||% 10000L))
    sec("use-community-specific-weighting",
        fi(cfg$community_specific_weighting %||% 0L))
    sec("number-of-base-mcmc-samples",
        fi(cfg$n_base_sampling %||% 500L))
    sec("number-of-burnin-mcmc-samples",
        fi(cfg$n_burnin_sampling %||% 10L))
    sec("number-of-burnin-mcmc-iterations",
        fi(cfg$n_burnin_iter %||% 10L))
    sec("number-of-samplings-for-mc-error",
        fi(cfg$n_sampling_for_mce %||% 0L))
    sec("use-bootstrap-for-mc-error",
        fi(cfg$use_bootstrap_for_mce %||% 0L))
    sec("do-not-calculate-average-variance-for-mc-error",
        fi(cfg$skip_Evar_for_mce %||% 1L))
    sec("check-missingness",  fi(cfg$check_missingness %||% 0L))
    sec("check-mixing",       fi(cfg$check_mixing      %||% 0L))
    sec("check-runtime",      fi(cfg$check_runtime     %||% 0L))

    # ---- 12.  Other scalar parameters --------------------------------------
    sec("relative-infectivity-of-asymptomatic-case-for-simulation",
        fd(cfg$asym_effect_sim %||% 0.5))
    sec("relative-infectivity-of-asymptomatic-case-for-estimation",
        fd(cfg$asym_effect_est %||% 0.5))
    sec("members-share-common-contact-history-within-communities",
        fi(cfg$common_contact_history_within_community %||% 1L))
    sec("automatically-generate-c2p-contact-file",
        fi(cfg$generate_c2p_contact %||% 1L))
    sec("automatically-generate-p2p-contact-file",
        fi(cfg$generate_p2p_contact %||% 1L))
    sec("use-c2p-offset",           fi(cfg$c2p_offset               %||% 0L))
    sec("use-p2p-offset",           fi(cfg$p2p_offset               %||% 0L))
    sec("adjust-for-selection-bias",fi(cfg$adjust_for_left_truncation %||% 1L))
    sec("adjust-for-right-censoring",fi(cfg$adjust_for_right_censoring %||% 1L))
    sec("prefix-index-cases",       fi(cfg$preset_index             %||% 1L))

    {
        use_b <- as.integer(cfg$use_index_cases_to_improve_b %||% 0L)
        if (isTRUE(use_b == 1L))
            sec("use-onset-days-of-index-cases-to-improve-estimation-of-b",
                paste0(fi(use_b), ":", fi(cfg$n_c2p_group %||% 0L)))
        else
            sec("use-onset-days-of-index-cases-to-improve-estimation-of-b",
                "0:")
    }

    sec("epidemic-duration-for-calculating-CPI",
        fi(cfg$CPI_duration %||% 0L))

    {
        prov <- as.integer(cfg$effective_bounds_provided %||% 0L)
        n_p  <- as.integer(cfg$n_p_mode %||% 0L)
        if (isTRUE(prov == 1L) && n_p > 0L) {
            body <- paste0(fi(prov), ":")
            for (i in seq_len(n_p))
                body <- c(body,
                          paste(fi(cfg$effective_lower_infectious[i]),
                                fi(cfg$effective_upper_infectious[i])))
            sec("effective-lower-upper-bounds-of-infectious-days-relative-to-symptom-onset-day",
                body)
        } else {
            sec("effective-lower-upper-bounds-of-infectious-days-relative-to-symptom-onset-day",
                "0:")
        }
    }

    # ---- 13.  SAR covariate sets -------------------------------------------
    {
        n_sets <- as.integer(cfg$SAR_n_covariate_sets %||% 0L)
        n_tic  <- as.integer(cfg$n_time_ind_covariate %||% 0L)
        n_tdc  <- as.integer(cfg$n_time_dep_covariate %||% 0L)
        n_p    <- as.integer(cfg$n_p_mode             %||% 0L)

        if (n_sets > 0L && n_p > 0L) {
            body <- paste0(fi(n_sets), ":")

            # TIC blocks — one sus-ind / inf-ind per set
            if (n_tic > 0L) {
                sti <- cfg$SAR_sus_time_ind_covariate
                iti <- cfg$SAR_inf_time_ind_covariate
                for (nn in seq_len(n_sets)) {
                    body <- c(body, "",
                              paste("sus-ind", fvd(sti[nn, ])),
                              paste("inf-ind", fvd(iti[nn, ])))
                }
            }

            # TDC blocks
            if (n_tdc > 0L) {
                lo  <- as.integer(cfg$SAR_time_dep_lower %||% 0L)
                hi  <- as.integer(cfg$SAR_time_dep_upper %||% 0L)
                len <- hi - lo + 1L
                body <- c(body, "", paste(fi(lo), fi(hi)))
                std  <- cfg$SAR_sus_time_dep_covariate
                itd  <- cfg$SAR_inf_time_dep_covariate

                for (nn in seq_len(n_sets)) {
                    body <- c(body, "")
                    # sus-dep: skip rows when all zero
                    if (is.null(std[[nn]]) || all(std[[nn]] == 0)) {
                        body <- c(body, "sus-dep 0")
                    } else {
                        body <- c(body, paste0("sus-dep ", fi(len)))
                        for (t in seq_len(len))
                            body <- c(body,
                                      paste(fi(lo + t - 1L), fi(lo + t - 1L),
                                            fvd(std[[nn]][t, ])))
                    }
                    # inf-dep
                    if (is.null(itd[[nn]]) || all(itd[[nn]] == 0)) {
                        body <- c(body, "inf-dep 0")
                    } else {
                        body <- c(body, paste0("inf-dep ", fi(len)))
                        for (t in seq_len(len))
                            body <- c(body,
                                      paste(fi(lo + t - 1L), fi(lo + t - 1L),
                                            fvd(itd[[nn]][t, ])))
                    }
                }
            }

            sec("covariates-for-calculating-SAR-provided", body)
        } else {
            sec("covariates-for-calculating-SAR-provided", "0:")
        }
    }

    # ---- 14.  R0 multiplier ------------------------------------------------
    {
        prov <- as.integer(cfg$R0_multiplier_provided %||% 0L)
        n_p  <- as.integer(cfg$n_p_mode %||% 0L)
        if (isTRUE(prov > 0L) && n_p > 0L) {
            body <- paste0(fi(prov), ":")
            for (i in seq_len(n_p))
                body <- c(body,
                          paste(fd(cfg$R0_multiplier[i]),
                                fd(cfg$R0_multiplier_var[i])))
            sec("multiplier-for-calculating-R0", body)
        } else {
            sec("multiplier-for-calculating-R0", "0:")
        }
    }

    {
        div <- as.integer(cfg$R0_divide_by_time %||% 0L)
        n_p <- as.integer(cfg$n_p_mode %||% 0L)
        if (isTRUE(div == 1L) && n_p > 0L) {
            body <- c(paste0(fi(div), ":"),
                      paste(fi(cfg$R0_divide_start), fi(cfg$R0_divide_stop)))
            if (n_p != 2L) body <- c(body, fi(cfg$R0_window_size))
            sec("serial-division-of-epidemic-for-calculating-time-varying-R0", body)
        } else {
            sec("serial-division-of-epidemic-for-calculating-time-varying-R0", "0:")
        }
    }

    # ---- 15.  Output / diagnostic switches ---------------------------------
    sec("goodness-of-fit",          fi(cfg$goodness_of_fit     %||% 0L))
    sec("perform-statistical-test", fi(cfg$stat_test           %||% 0L))
    sec("do-not-estimate-variance", fi(cfg$skip_variance       %||% 0L))
    sec("print-covariance-matrix",  fi(cfg$print_covariance    %||% 0L))
    sec("do-not-output-estimates",  fi(cfg$skip_output         %||% 0L))
    sec("simplify-output",          fi(cfg$simplify_output     %||% 0L))
    sec("simplify-output-SAR",      fi(cfg$simplify_output_SAR %||% 0L))
    sec("simplify-output-R0",       fi(cfg$simplify_output_R0  %||% 0L))
    sec("run-transtat-silently",    fi(cfg$silent_run          %||% 1L))
    sec("write-error-log",          fi(cfg$write_error_log     %||% 0L))

    # ---- Write file --------------------------------------------------------
    writeLines(lns, out_path)
    invisible(out_path)
}
