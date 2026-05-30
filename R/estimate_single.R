#' Maximum Likelihood Estimation for a Single Epidemic Dataset
#'
#' Runs the TranStat chain-binomial MLE engine on observed or simulated
#' epidemic data.  Population data frames are passed directly to C without
#' intermediate file I/O, making this the preferred low-level estimation
#' routine for simulation studies and bootstrapping.
#'
#' @section Parameter scales:
#' Transmission probabilities (\code{b}, \code{p}) and group parameters
#' (\code{u}, \code{q}) are estimated on the logit scale and
#' back-transformed before reporting.  Covariate coefficients are estimated
#' on the log scale and reported as odds ratios.  Confidence intervals are
#' obtained by back-transforming the Wald interval on the logit / log scale.
#' The Wald z-statistic tests \eqn{\mathrm{logit}(p) = 0} (i.e., \eqn{p =
#' 0.5}) for probability parameters and \eqn{\log(\mathrm{OR}) = 0}
#' (i.e., \eqn{\mathrm{OR} = 1}) for covariate coefficients.
#'
#' @param data_list Named list as returned by \code{\link{gen_population}},
#'   \code{\link{simulate_single}}, or \code{\link{read_population}}.  Must
#'   contain \code{$pop}, \code{$community}, \code{$time_ind_covariate}, and
#'   \code{$time_dep_covariate}.  Optional elements \code{$c2p_contact},
#'   \code{$p2p_contact}, and \code{$impute} are used when non-\code{NULL}.
#' @param cfg Named list returned by \code{\link{read_config}}.
#' @param i_inc Integer (1-based). Incubation-period group index. Default
#'   \code{1L}.
#' @param i_inf Integer (1-based). Infectious-period group index. Default
#'   \code{1L}.
#' @param seed Integer. RNG seed for MCMC / MCEM state-space sampling.
#'   Default \code{12345678L}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{estimates}}{Data frame with one row per model parameter
#'     (columns: \code{parameter}, \code{estimate}, \code{se},
#'     \code{ci_lower}, \code{ci_upper}, \code{z}, \code{p_value}).}
#'   \item{\code{log_likelihood}}{Log-likelihood at the MLE.}
#'   \item{\code{error_code}}{Integer returned by the C estimator
#'     (0 = converged successfully).}
#' }
#'
#' @seealso \code{\link{gen_population}}, \code{\link{simulate_single}},
#'   \code{\link{read_config}}, \code{\link{gen_impute}}
#'
#' @examples
#' \dontrun{
#' cfg <- read_config(system.file("extdata", "CaseStudy1", "config.file",
#'                                package = "ChainBinomial"))
#' pop <- gen_population(n_community = 200, community_size = 5,
#'                       day_epi_stop = 14, case_ascertained = 1L)
#' sim <- simulate_single(dat, cfg, seed = 42L)
#' fit <- estimate_single(sim, cfg, seed = 42L)
#' fit$estimates
#' fit$log_likelihood
#' }
#'
#' @export
estimate_single <- function(data_list, cfg,
                             i_inc = 1L, i_inf = 1L,
                             seed  = 12345678L) {

    stopifnot(is.list(data_list), is.list(cfg))
    stopifnot(i_inc >= 1L, i_inc <= cfg$n_inc)
    stopifnot(i_inf >= 1L, i_inf <= cfg$n_inf)

    cfg2    <- cfg
    n_peq   <- cfg$n_par_equiclass

    # -----------------------------------------------------------------------
    # Identify which parameter positions use logit (prob) vs log (coeff).
    # Positions 0 .. n_prob_pars-1 → logit.
    # -----------------------------------------------------------------------
    n_prob_pars <- cfg$n_b_mode + cfg$n_p_mode + cfg$n_u_mode + cfg$n_q_mode

    # Transform one raw value given its 0-based parameter position.
    tr_one <- function(pos0, raw) {
        if (pos0 < n_prob_pars) {
            v <- max(1e-10, min(1 - 1e-10, raw))   # clamp before logit
            log(v / (1.0 - v))
        } else {
            v <- max(1e-300, raw)                   # clamp before log
            log(v)
        }
    }

    # -----------------------------------------------------------------------
    # Transform ini_par_effective (n_ini × n_peq, raw → logit/log).
    # Column j corresponds to equivalence class j; the 0-based position of
    # its first member determines which transform applies.
    # Flatten row-major for C.
    # -----------------------------------------------------------------------
    ini_tr <- cfg$ini_par_effective   # already a matrix (n_ini × n_peq)
    if (!is.null(cfg$par_equiclass) && n_peq > 0L &&
            isTRUE(cfg$ini_par_provided == 1L)) {
        for (j in seq_len(n_peq)) {
            pos0 <- cfg$par_equiclass[[j]]$member[1L] - 1L
            ini_tr[, j] <- vapply(ini_tr[, j], tr_one, double(1), pos0 = pos0)
        }
    }
    # When not provided, leave as 0.0 (= logit(0.5) / log(1) — neutral start)
    cfg2$ini_par_effective_flat <- as.vector(t(ini_tr))   # row-major

    # -----------------------------------------------------------------------
    # Transform search bounds (one per equivalence class, raw → logit/log).
    # When not provided the values default to 0.0; C ignores them in that case.
    # -----------------------------------------------------------------------
    lower_tr <- cfg$lower_search_bound
    upper_tr <- cfg$upper_search_bound
    if (!is.null(cfg$par_equiclass) && n_peq > 0L &&
            isTRUE(cfg$search_bound_provided == 1L)) {
        for (j in seq_len(n_peq)) {
            pos0 <- cfg$par_equiclass[[j]]$member[1L] - 1L
            lower_tr[j] <- tr_one(pos0, lower_tr[j])
            upper_tr[j] <- tr_one(pos0, upper_tr[j])
        }
    }
    cfg2$lower_search_bound_tr <- lower_tr
    cfg2$upper_search_bound_tr <- upper_tr

    # -----------------------------------------------------------------------
    # Transform fixed parameter values (indexed by their raw parameter position).
    # -----------------------------------------------------------------------
    fpv_tr <- numeric(0)
    if (cfg$n_par_fixed > 0L) {
        fpv_tr <- numeric(cfg$n_par_fixed)
        for (k in seq_len(cfg$n_par_fixed)) {
            pos0      <- cfg$par_fixed_id[k] - 1L   # 0-based
            fpv_tr[k] <- tr_one(pos0, cfg$par_fixed_value[k])
        }
    }
    cfg2$par_fixed_value_tr <- fpv_tr

    # -----------------------------------------------------------------------
    # Effective infectious bounds used for SAR / R0 computation.
    # Default: use the limits of the selected infectious-period group.
    # -----------------------------------------------------------------------
    if (is.null(cfg$effective_lower_infectious) || cfg$n_p_mode == 0L) {
        cfg2$effective_lower_infectious_tr <-
            rep(as.integer(cfg$lower_inf[i_inf]), max(cfg$n_p_mode, 1L))
        cfg2$effective_upper_infectious_tr <-
            rep(as.integer(cfg$upper_inf[i_inf]), max(cfg$n_p_mode, 1L))
    } else {
        cfg2$effective_lower_infectious_tr <-
            as.integer(cfg$effective_lower_infectious)
        cfg2$effective_upper_infectious_tr <-
            as.integer(cfg$effective_upper_infectious)
    }

    # -----------------------------------------------------------------------
    # NA → 0 for integer fields the C layer reads unconditionally.
    # -----------------------------------------------------------------------
    if (is.na(cfg2$PreIllness_covariate_id)) cfg2$PreIllness_covariate_id <- 0L
    if (is.na(cfg2$Illness_covariate_id))    cfg2$Illness_covariate_id    <- 0L
    if (is.na(cfg2$RxIllness_covariate_id))  cfg2$RxIllness_covariate_id  <- 0L
    if (is.na(cfg2$RxIllness_prob))          cfg2$RxIllness_prob          <- 0.0
    if (is.na(cfg2$RxIllness_duration))      cfg2$RxIllness_duration      <- 0L
    if (is.na(cfg2$RxIllness_index_only))    cfg2$RxIllness_index_only    <- 0L

    # -----------------------------------------------------------------------
    # NULL → empty numeric(0) for optional R0 fields.
    # -----------------------------------------------------------------------
    if (is.null(cfg2$R0_multiplier))     cfg2$R0_multiplier     <- numeric(0)
    if (is.null(cfg2$R0_multiplier_var)) cfg2$R0_multiplier_var <- numeric(0)

    # -----------------------------------------------------------------------
    # Coerce columns that the C layer reads with REAL() to double.
    # read.table infers integer for whole-number columns; C requires double.
    #   pop$weight        → col 12 (1-based), read with REAL()
    #   TIC value cols    → cols 2..n_tic+1, read with REAL()
    #   TDC value cols    → cols 4..n_tdc+3, read with REAL()
    # -----------------------------------------------------------------------
    pop_df <- data_list$pop
    pop_df[["weight"]] <- as.double(pop_df[["weight"]])

    tic <- data_list$time_ind_covariate
    if (!is.null(tic) && cfg$n_time_ind_covariate > 0L) {
        for (j in seq_len(cfg$n_time_ind_covariate)) {
            tic[[j + 1L]] <- as.double(tic[[j + 1L]])
        }
    }
    tdc <- data_list$time_dep_covariate
    if (!is.null(tdc) && cfg$n_time_dep_covariate > 0L) {
        for (j in seq_len(cfg$n_time_dep_covariate)) {
            tdc[[j + 3L]] <- as.double(tdc[[j + 3L]])
        }
    }

    # C reads the contact offset column with REAL().
    # Shared history: community_id start_day stop_day contact_mode offset ignore
    #                 col 5 (1-based) = offset
    # Individualised p2p: start_day stop_day person_i person_j contact_mode offset ignore
    #                     col 6 (1-based) = offset
    c2p <- data_list$c2p_contact
    if (!is.null(c2p)) {
        c2p[[5L]] <- as.double(c2p[[5L]])   # offset col (both shared and individualised)
    }
    p2p <- data_list$p2p_contact
    if (!is.null(p2p)) {
        off_col <- if (cfg$common_contact_history_within_community == 1L) 5L else 6L
        p2p[[off_col]] <- as.double(p2p[[off_col]])
    }

    # -----------------------------------------------------------------------
    # Call the C estimation engine.
    # NULL list elements (c2p_contact, p2p_contact, impute) become
    # R_NilValue in C, which the C function handles gracefully.
    # -----------------------------------------------------------------------
    raw <- .Call("r_estimate_single",
                 pop_df,
                 tic,
                 data_list$community,
                 tdc,
                 c2p,
                 p2p,
                 data_list$impute,
                 cfg2,
                 as.integer(i_inc),
                 as.integer(i_inf),
                 as.integer(seed))

    # -----------------------------------------------------------------------
    # Post-processing.
    # raw$est      : double[n_par]   already on raw (probability / OR) scale
    # raw$var      : double[n_par²]  col-major covariance, raw scale
    # raw$var_logit: double[n_par²]  col-major covariance, logit / log scale
    # -----------------------------------------------------------------------
    n_par   <- cfg$n_par
    est_raw <- raw$est                                              # raw scale from C
    V_raw   <- matrix(raw$var,       nrow = n_par, ncol = n_par)   # col-major → matrix
    V_tr    <- matrix(raw$var_logit, nrow = n_par, ncol = n_par)   # col-major → matrix

    # Per-parameter logical: TRUE = probability param (logit scale in var_logit)
    is_prob <- (seq_len(n_par) - 1L) < n_prob_pars

    # Logit / log of the MLE (needed for Wald CI and z-statistic).
    # Clamp before transforming to avoid log(0) or logit(0/1).
    est_ll <- numeric(n_par)
    for (k in seq_len(n_par)) {
        if (is_prob[k]) {
            v         <- max(1e-300, min(1 - 1e-300, est_raw[k]))
            est_ll[k] <- log(v / (1.0 - v))                # logit
        } else {
            est_ll[k] <- log(max(1e-300, est_raw[k]))       # log
        }
    }

    # SE on raw scale from the raw-scale variance matrix
    se_raw <- sqrt(pmax(0.0, diag(V_raw)))

    # SE on logit / log scale from var_logit diagonal
    se_tr  <- sqrt(pmax(0.0, diag(V_tr)))

    # 95% CI: Wald interval on logit/log scale, back-transformed to raw scale
    z95   <- qnorm(0.975)
    lo_ll <- est_ll - z95 * se_tr
    hi_ll <- est_ll + z95 * se_tr
    ci_lo <- ifelse(is_prob, 1.0 / (1.0 + exp(-lo_ll)), exp(lo_ll))
    ci_hi <- ifelse(is_prob, 1.0 / (1.0 + exp(-hi_ll)), exp(hi_ll))

    # Wald z-statistic (logit/log scale) — internal only, not in output
    z_stat  <- ifelse(se_tr > 0.0, est_ll / se_tr, NA_real_)
    # Two-sided p-value; set to NA for base transmission parameters below
    p_value <- 2.0 * pnorm(-abs(z_stat))

    # -----------------------------------------------------------------------
    # Parameter names.
    #
    # Base parameters (b, p, u, q) use the "option-2" convention that matches
    # cfg$par_labels: the trailing digit is omitted when there is only one mode
    # (e.g. "b" not "b1", but "b1","b2" when n_b_mode > 1).
    #
    # Covariate parameters are labelled using the actual column names carried in
    # data_list$time_ind_covariate / data_list$time_dep_covariate, so that
    # user-supplied names (e.g. "age", "antiviral") flow through to all output.
    # The fallback when those data frames are unavailable is "x{global_index}".
    # -----------------------------------------------------------------------

    # Helper: option-2 base-parameter label vector
    .lbl <- function(prefix, n) {
        if (n == 0L) character(0L)
        else if (n == 1L) prefix
        else paste0(prefix, seq_len(n))
    }

    # Helper: global covariate index k (1-based) → column name from data_list
    .n_tic <- cfg$n_time_ind_covariate
    cov_label <- function(k) {
        if (k <= .n_tic) {
            tic <- data_list$time_ind_covariate
            if (!is.null(tic) && ncol(tic) >= k + 1L)
                return(names(tic)[k + 1L])
        } else {
            tdc <- data_list$time_dep_covariate
            col <- k - .n_tic + 3L          # cols 1-3 are id/day_start/day_stop
            if (!is.null(tdc) && ncol(tdc) >= col)
                return(names(tdc)[col])
        }
        paste0("x", k)                      # fallback
    }

    n_base <- cfg$n_b_mode + cfg$n_p_mode + cfg$n_u_mode + cfg$n_q_mode

    # p_value is only meaningful for covariate effects (testing H0: OR = 1).
    # Base transmission parameters (b, p, u, q) get NA.
    if (n_base > 0L) p_value[seq_len(n_base)] <- NA_real_

    par_names <- c(
        .lbl("b", cfg$n_b_mode),
        .lbl("p", cfg$n_p_mode),
        .lbl("u", cfg$n_u_mode),
        .lbl("q", cfg$n_q_mode),
        if (cfg$n_c2p_covariate      > 0L)
            paste0("c2p_",
                   vapply(cfg$c2p_covariate[seq_len(cfg$n_c2p_covariate)],
                          cov_label, character(1L))),
        if (cfg$n_sus_p2p_covariate  > 0L)
            paste0("p2p_s_",
                   vapply(cfg$sus_p2p_covariate[seq_len(cfg$n_sus_p2p_covariate)],
                          cov_label, character(1L))),
        if (cfg$n_inf_p2p_covariate  > 0L)
            paste0("p2p_i_",
                   vapply(cfg$inf_p2p_covariate[seq_len(cfg$n_inf_p2p_covariate)],
                          cov_label, character(1L))),
        if (!is.null(cfg$interaction) && cfg$n_int_p2p_covariate > 0L)
            vapply(seq_len(cfg$n_int_p2p_covariate), function(j)
                paste0("p2p_int_",
                       cov_label(cfg$interaction[[j]][1L]),
                       "_",
                       cov_label(cfg$interaction[[j]][2L])), character(1L)),
        .lbl("pat", cfg$n_pat_covariate),
        .lbl("imm", cfg$n_imm_covariate)
    )

    estimates <- data.frame(
        parameter = par_names,
        estimate  = est_raw,
        se        = se_raw,
        ci_lower  = ci_lo,
        ci_upper  = ci_hi,
        p_value   = p_value,
        row.names = NULL,
        stringsAsFactors = FALSE
    )

    # Attach parameter names to rows/cols of both matrices for readability
    dimnames(V_raw) <- list(par_names, par_names)
    dimnames(V_tr)  <- list(par_names, par_names)

    # -----------------------------------------------------------------------
    # Unadjusted secondary attack rates (SAR) and basic reproduction number
    # (R0).  Mirrors core.h lines 1545-1699.
    #
    # For each p-to-p mode k the SAR is the probability that a susceptible
    # contact is infected by a single infective over the effective infectious
    # period:
    #
    #   SAR_k = 1 - prod_{t=eff_lo_k}^{eff_hi_k} (1 - p_k * s_t)
    #
    # where s_t = prob_infectious[t - lower_inf + 1] is the daily
    # infectiousness weight.  The delta-method variance uses the full
    # var_logit matrix so that covariate parameters (which influence the
    # p2p probability in the adjusted SAR) can contribute in future
    # extensions; here only the p_k parameter has a non-zero derivative.
    #
    # Confidence intervals are Wald intervals on the logit(SAR) scale,
    # back-transformed.
    #
    # R0 = sum_k  R0_multiplier_k * SAR_k   (when R0_multiplier_provided).
    # Its CI uses the log-normal approximation: exp(log(R0) ± 1.96*SE/R0).
    # An additional variance term accounts for uncertainty in the multiplier
    # itself (R0_multiplier_var).
    #
    # Note on R0 SE vs. the reference transtat output:
    #   core.h computes unadjusted SAR derivatives into der_mat, then
    #   (when SAR_n_covariate_sets > 0) overwrites der_mat with the
    #   covariate-adjusted SAR derivatives before the unadjusted R0
    #   variance loop.  Because r_estimate_single sets
    #   SAR_n_covariate_sets = 0, that overwrite never occurs here, and
    #   the R0 variance is computed from the unadjusted SAR derivatives —
    #   which is the mathematically correct formula for the unadjusted R0.
    # -----------------------------------------------------------------------
    n_b    <- cfg$n_b_mode
    n_p    <- cfg$n_p_mode
    sar_df <- NULL
    r0_df  <- NULL

    if (n_p > 0L) {
        prob_inf  <- cfg$prob_inf[[i_inf]]        # daily infectiousness weights
        lower_inf <- cfg$lower_inf[i_inf]          # offset into prob_inf
        eff_lo    <- cfg2$effective_lower_infectious_tr   # length n_p
        eff_hi    <- cfg2$effective_upper_infectious_tr

        SAR0        <- numeric(n_p)
        se_SAR0     <- numeric(n_p)
        lo_SAR0     <- numeric(n_p)
        hi_SAR0     <- numeric(n_p)
        # der_mat_sar[k, j] = d(SAR0_k) / d(logit-param_j)
        der_mat_sar <- matrix(0.0, nrow = n_p, ncol = n_par)

        for (k in seq_len(n_p)) {
            p_esc <- 1.0
            der   <- numeric(n_par)     # d(log p_esc) / d(logit-param_j)
            lp_k  <- est_ll[n_b + k]   # logit of p[k]; 1-based: position n_b+k
            p_k   <- 1.0 / (1.0 + exp(-lp_k))  # inv_logit(lp_k)

            t_seq <- seq.int(eff_lo[k], eff_hi[k])
            for (t in t_seq) {
                l  <- t - lower_inf + 1L          # 1-based index into prob_inf
                s  <- prob_inf[l]
                ff <- p_k * s                      # = inv_logit(lp_k) * s
                f  <- 1.0 - ff
                # d(log f)/d(lp_k) = -ff*(1 - p_k)/f  (= -ff*(1 - ff/s)/f)
                log_f_lp_k        <- -ff * (1.0 - p_k) / f
                der[n_b + k]      <- der[n_b + k] + log_f_lp_k
                p_esc             <- p_esc * f
            }

            SAR0[k]          <- 1.0 - p_esc
            der_mat_sar[k, ] <- -p_esc * der       # d(SAR0_k)/d(logit-param_j)

            # Variance via delta method: d' * var_logit * d
            d       <- der_mat_sar[k, ]
            var_sar <- as.numeric(d %*% V_tr %*% d)
            se_SAR0[k] <- sqrt(max(0.0, var_sar))

            # 95% CI on logit(SAR) scale, then back-transform
            sar_k  <- SAR0[k]
            sar_k  <- max(1e-300, min(1 - 1e-300, sar_k))   # guard logit
            se_logit_sar <- se_SAR0[k] / (sar_k * (1.0 - sar_k))
            logit_sar    <- log(sar_k / (1.0 - sar_k))
            lo_SAR0[k]   <- 1.0 / (1.0 + exp(-(logit_sar - 1.96 * se_logit_sar)))
            hi_SAR0[k]   <- 1.0 / (1.0 + exp(-(logit_sar + 1.96 * se_logit_sar)))
        }

        group_labels <- if (n_p == 1L) "overall" else paste0("p", seq_len(n_p))
        sar_df <- data.frame(
            group    = group_labels,
            SAR      = SAR0,
            se       = se_SAR0,
            ci_lower = lo_SAR0,
            ci_upper = hi_SAR0,
            row.names = NULL,
            stringsAsFactors = FALSE
        )

        # ---- R0 multiplier (defined here so adjusted R0 block can reuse) ----
        mult   <- NULL
        mult_v <- NULL
        if (isTRUE(cfg$R0_multiplier_provided > 0L) &&
                !is.null(cfg$R0_multiplier) &&
                length(cfg$R0_multiplier) >= n_p) {
            mult   <- cfg$R0_multiplier[seq_len(n_p)]
            mult_v <- cfg$R0_multiplier_var[seq_len(n_p)]
        }

        # ---- Unadjusted R0 ----
        if (!is.null(mult)) {
            R0_val <- sum(mult * SAR0)

            # der1[j] = d(R0)/d(logit-param_j) = sum_k mult_k * der_mat_sar[k,j]
            der1 <- as.numeric(mult %*% der_mat_sar)

            # Variance: quadratic form in var_logit + multiplier uncertainty
            var_R0 <- as.numeric(der1 %*% V_tr %*% der1) + sum(SAR0^2 * mult_v)
            se_R0  <- sqrt(max(0.0, var_R0))

            # CI on log(R0) scale (log-normal), back-transformed
            lo_R0 <- exp(log(R0_val) - 1.96 * se_R0 / R0_val)
            hi_R0 <- exp(log(R0_val) + 1.96 * se_R0 / R0_val)

            r0_df <- data.frame(
                R0       = R0_val,
                se       = se_R0,
                ci_lower = lo_R0,
                ci_upper = hi_R0,
                row.names = NULL
            )
        }

        # -----------------------------------------------------------------------
        # Covariate-adjusted SAR and R0 (mirrors core.h lines 1583-1638 and
        # 1674-1698).
        #
        # For each covariate set m and p-mode k the adjusted SAR uses
        #   logit_f = logit(p_k) + sum_j p2p_cov[j] * log_OR[j]
        # where p2p_cov is built by organize_p2p_covariate_4SAR() in C,
        # replicated here in .org_p2p_cov().
        #
        # The delta-method gradient (der_mat_adj) is re-computed for every m.
        # After the loop it holds the LAST set's derivatives; C reuses these
        # for all adjusted R0 SEs (matching the reference output).
        # -----------------------------------------------------------------------

        # --- R equivalent of organize_p2p_covariate_4SAR ---
        .org_p2p_cov <- function(m0, r0) {
            # m0: 0-based set index; r0: 0-based time index
            n_tic <- cfg$n_time_ind_covariate
            n_tdc <- cfg$n_time_dep_covariate
            n_cov <- n_tic + n_tdc
            sus_v <- numeric(n_cov)
            inf_v <- numeric(n_cov)
            if (n_tic > 0L && !is.null(cfg$SAR_sus_time_ind_covariate)) {
                sus_v[seq_len(n_tic)] <- cfg$SAR_sus_time_ind_covariate[m0 + 1L, ]
                inf_v[seq_len(n_tic)] <- cfg$SAR_inf_time_ind_covariate[m0 + 1L, ]
            }
            if (n_tdc > 0L && !is.null(cfg$SAR_inf_time_dep_covariate)) {
                r1  <- r0 + 1L
                len <- nrow(cfg$SAR_inf_time_dep_covariate[[m0 + 1L]])
                if (r1 >= 1L && r1 <= len) {
                    sus_v[n_tic + seq_len(n_tdc)] <-
                        cfg$SAR_sus_time_dep_covariate[[m0 + 1L]][r1, ]
                    inf_v[n_tic + seq_len(n_tdc)] <-
                        cfg$SAR_inf_time_dep_covariate[[m0 + 1L]][r1, ]
                }
            }
            p2p_v <- numeric(cfg$n_p2p_covariate)
            l     <- 0L
            if (cfg$n_sus_p2p_covariate > 0L)
                for (j in seq_len(cfg$n_sus_p2p_covariate)) {
                    l <- l + 1L; p2p_v[l] <- sus_v[cfg$sus_p2p_covariate[j]]
                }
            if (cfg$n_inf_p2p_covariate > 0L)
                for (j in seq_len(cfg$n_inf_p2p_covariate)) {
                    l <- l + 1L; p2p_v[l] <- inf_v[cfg$inf_p2p_covariate[j]]
                }
            if (!is.null(cfg$interaction) && cfg$n_int_p2p_covariate > 0L)
                for (j in seq_len(cfg$n_int_p2p_covariate)) {
                    l <- l + 1L
                    p2p_v[l] <- sus_v[cfg$interaction[[j]][1L]] *
                                 inf_v[cfg$interaction[[j]][2L]]
                }
            p2p_v
        }

        sar_adj_df <- NULL
        r0_adj_df  <- NULL

        has_adj <- isTRUE(cfg$SAR_n_covariate_sets > 0L) &&
                   (!is.null(cfg$SAR_inf_time_dep_covariate) ||
                    !is.null(cfg$SAR_inf_time_ind_covariate))

        if (has_adj) {
            n_sets  <- cfg$SAR_n_covariate_sets
            n_p2p   <- cfg$n_p2p_covariate

            # Positions (1-based) of p2p parameters in the full parameter vector
            p2p_off <- n_b + n_p + cfg$n_u_mode + cfg$n_q_mode + cfg$n_c2p_covariate
            p2p_idx <- p2p_off + seq_len(n_p2p)
            # log-OR values (= logit-scale coefficients for p2p parameters)
            coeff_p2p <- if (n_p2p > 0L) est_ll[p2p_idx] else numeric(0L)

            SAR_adj     <- matrix(0.0, nrow = n_sets, ncol = n_p)
            se_SAR_adj  <- matrix(0.0, nrow = n_sets, ncol = n_p)
            lo_SAR_adj  <- matrix(0.0, nrow = n_sets, ncol = n_p)
            hi_SAR_adj  <- matrix(0.0, nrow = n_sets, ncol = n_p)
            der_mat_adj <- matrix(0.0, nrow = n_p, ncol = n_par)  # overwritten per m

            for (m in seq_len(n_sets)) {       # m: 1-based in R, 0-based in C
                for (k in seq_len(n_p)) {
                    p_esc <- 1.0
                    der   <- numeric(n_par)
                    lp_k  <- est_ll[n_b + k]

                    for (t in seq.int(eff_lo[k], eff_hi[k])) {
                        r0     <- t - cfg$SAR_time_dep_lower   # 0-based time index
                        p2p_cv <- .org_p2p_cov(m - 1L, r0)

                        cov_eff <- if (n_p2p > 0L) sum(p2p_cv * coeff_p2p) else 0.0
                        lf      <- lp_k + cov_eff
                        l_idx   <- t - lower_inf + 1L
                        s_t     <- prob_inf[l_idx]
                        inv_lf  <- 1.0 / (1.0 + exp(-lf))
                        ff      <- inv_lf * s_t
                        f       <- 1.0 - ff
                        lf_lpk  <- -ff * (1.0 - inv_lf) / f  # d(log f)/d(lp_k)

                        der[n_b + k] <- der[n_b + k] + lf_lpk
                        if (n_p2p > 0L)
                            der[p2p_idx] <- der[p2p_idx] + lf_lpk * p2p_cv
                        p_esc <- p_esc * f
                    }

                    SAR_adj[m, k]    <- 1.0 - p_esc
                    der_mat_adj[k, ] <- -p_esc * der

                    d_k     <- der_mat_adj[k, ]
                    var_sar <- as.numeric(d_k %*% V_tr %*% d_k)
                    se_SAR_adj[m, k] <- sqrt(max(0.0, var_sar))

                    sar_k    <- max(1e-300, min(1.0 - 1e-300, SAR_adj[m, k]))
                    se_lgit  <- se_SAR_adj[m, k] / (sar_k * (1.0 - sar_k))
                    lgit_sar <- log(sar_k / (1.0 - sar_k))
                    lo_SAR_adj[m, k] <- 1.0 / (1.0 + exp(-(lgit_sar - 1.96 * se_lgit)))
                    hi_SAR_adj[m, k] <- 1.0 / (1.0 + exp(-(lgit_sar + 1.96 * se_lgit)))
                }
            }
            # der_mat_adj now holds derivatives from the LAST covariate set.
            # C reuses this same der_mat for all adjusted R0 SEs (see core.h
            # lines 1682-1683 — der_mat is not reset between R0_adj iterations).

            sar_adj_rows <- lapply(seq_len(n_sets), function(m) {
                data.frame(
                    covariate_set = m - 1L,
                    group         = group_labels,
                    SAR           = SAR_adj[m, ],
                    se            = se_SAR_adj[m, ],
                    ci_lower      = lo_SAR_adj[m, ],
                    ci_upper      = hi_SAR_adj[m, ],
                    row.names     = NULL,
                    stringsAsFactors = FALSE
                )
            })
            sar_adj_df <- do.call(rbind, sar_adj_rows)

            # ---- Adjusted R0 ----
            if (!is.null(mult)) {
                # der1 uses der_mat_adj from LAST set (matches C behaviour)
                der1_last <- as.numeric(mult %*% der_mat_adj)

                R0_adj_val <- numeric(n_sets)
                se_R0_adj  <- numeric(n_sets)
                lo_R0_adj  <- numeric(n_sets)
                hi_R0_adj  <- numeric(n_sets)

                for (m in seq_len(n_sets)) {
                    R0_adj_val[m] <- sum(mult * SAR_adj[m, ])
                    der2          <- SAR_adj[m, ]
                    # var = quadratic form (last-set der1) + multiplier variance
                    var_R0_adj   <- as.numeric(der1_last %*% V_tr %*% der1_last) +
                                    sum(der2^2 * mult_v)
                    se_R0_adj[m] <- sqrt(max(0.0, var_R0_adj))
                    lo_R0_adj[m] <- exp(log(R0_adj_val[m]) -
                                        1.96 * se_R0_adj[m] / R0_adj_val[m])
                    hi_R0_adj[m] <- exp(log(R0_adj_val[m]) +
                                        1.96 * se_R0_adj[m] / R0_adj_val[m])
                }

                r0_adj_df <- data.frame(
                    covariate_set = seq_len(n_sets) - 1L,
                    R0            = R0_adj_val,
                    se            = se_R0_adj,
                    ci_lower      = lo_R0_adj,
                    ci_upper      = hi_R0_adj,
                    row.names     = NULL
                )
            }
        }
    }

    list(
        estimates      = estimates,
        SAR            = sar_df,                           # unadjusted SAR per p2p group
        SAR_adjusted   = sar_adj_df,                       # covariate-adjusted SAR (or NULL)
        R0             = r0_df,                            # unadjusted R0 (or NULL)
        R0_adjusted    = r0_adj_df,                        # covariate-adjusted R0 (or NULL)
        est_tr         = setNames(est_raw, par_names),     # direct C output (raw scale)
        var            = V_raw,                            # covariance on prob/OR scale
        var_logit      = V_tr,                             # covariance on logit/log scale
        log_likelihood = raw$log_likelihood,
        error_code     = raw$error_code
    )
}
