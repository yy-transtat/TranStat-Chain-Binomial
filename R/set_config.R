#' Set or Modify Configuration Parameters
#'
#' Creates a modified copy of an existing configuration object.  Only
#' parameters for which the caller supplies a non-\code{NULL} argument are
#' updated; everything else is left exactly as it was in \code{cfg}.
#'
#' @section Covariate specifications:
#' \code{c2p_covariate}, \code{sus_p2p_covariate}, \code{inf_p2p_covariate},
#' \code{pat_covariate}, and \code{imm_covariate} may be supplied as either
#' \itemize{
#'   \item an integer vector of global covariate indices (1-based, counting
#'     time-independent covariates first then time-dependent), or
#'   \item a character vector of covariate column names as they appear in
#'     \code{data_list$time_ind_covariate} or
#'     \code{data_list$time_dep_covariate}.  \code{data_list} must be
#'     supplied when names are used; an error is thrown for unrecognised
#'     names.
#' }
#'
#' @section Equivalence classes (\code{par_equiclass}):
#' Each element of the list must be a named list with a \code{member} field
#' that is either an integer vector of sequential parameter indices (1-based)
#' or a character vector of parameter labels (e.g.\
#' \code{"b"}, \code{"c2p_x1"}, \code{"p2p_i_antiviral"}).  When labels are
#' used \code{data_list} is recommended so that covariate labels resolve to
#' meaningful names.  The \code{size} field is always recomputed from
#' \code{length(member)}.
#'
#' @param cfg Named list returned by \code{\link{read_config}}.  Must not be
#'   \code{NULL}.
#' @param data_list Named list returned by \code{\link{read_population}} or
#'   \code{\link{gen_population}}.  When supplied, several counts are inferred
#'   automatically and covariate names can be used in place of integer indices.
#' @param path_in Character. Input data directory.
#' @param path_out Character. Output directory.
#' @param min_inc Integer vector (length \code{n_inc}).  Minimum incubation
#'   days per group.  Must be supplied together with \code{max_inc} and
#'   \code{prob_inc}.
#' @param max_inc Integer vector (length \code{n_inc}).  Maximum incubation
#'   days per group.
#' @param prob_inc List of numeric vectors (length \code{n_inc}), one
#'   probability vector per group.
#' @param lower_inf Integer vector (length \code{n_inf}).  Lower infectious
#'   bound per group.  Must be supplied together with \code{upper_inf} and
#'   \code{prob_inf}.
#' @param upper_inf Integer vector (length \code{n_inf}).  Upper infectious
#'   bound per group.
#' @param prob_inf List of numeric vectors (length \code{n_inf}), one
#'   probability vector per group.
#' @param c2p_covariate Integer or character vector.  Covariates affecting
#'   susceptibility for community-to-person transmission.
#' @param sus_p2p_covariate Integer or character vector.  Covariates affecting
#'   susceptibility for person-to-person transmission.
#' @param inf_p2p_covariate Integer or character vector.  Covariates affecting
#'   infectiousness for person-to-person transmission.
#' @param pat_covariate Integer or character vector.  Covariates affecting
#'   pathogenicity.
#' @param imm_covariate Integer or character vector.  Covariates affecting
#'   pre-season immunity.
#' @param par_equiclass List of equivalence classes.  Each element is a list
#'   with a \code{member} field (integer or character vector).
#' @param \dots Additional named arguments corresponding to any other element
#'   of \code{cfg}.  Values are assigned directly without validation.
#'
#' @return A modified copy of \code{cfg} with \code{par_labels} updated when
#'   \code{data_list} is provided.
#'
#' @seealso \code{\link{read_config}}, \code{\link{read_population}},
#'   \code{\link{update_par_labels}}, \code{\link{write_config}}
#'
#' @examples
#' \dontrun{
#' cfg       <- read_config(config_file)
#' data_list <- read_population(cfg, names_tdp = c("antiviral_A", "antiviral_B"))
#'
#' # Specify inf_p2p covariates by name
#' cfg2 <- set_config(cfg, data_list,
#'                    inf_p2p_covariate = c("antiviral_A", "antiviral_B"))
#'
#' # Equivalently, by global index
#' cfg2 <- set_config(cfg, data_list, inf_p2p_covariate = c(51L, 52L))
#'
#' # Equivalence class: c2p and p2p_s share the same value
#' cfg3 <- set_config(cfg2, data_list,
#'   par_equiclass = list(
#'     list(member = "b"),
#'     list(member = "p"),
#'     list(member = c("c2p_antiviral_A", "p2p_s_antiviral_A")),
#'     list(member = "p2p_i_antiviral_B")
#'   ))
#' }
#'
#' @export
set_config <- function(cfg,
                       data_list         = NULL,
                       path_in           = NULL,
                       path_out          = NULL,
                       min_inc           = NULL,
                       max_inc           = NULL,
                       prob_inc          = NULL,
                       lower_inf         = NULL,
                       upper_inf         = NULL,
                       prob_inf          = NULL,
                       c2p_covariate     = NULL,
                       sus_p2p_covariate = NULL,
                       inf_p2p_covariate = NULL,
                       pat_covariate     = NULL,
                       imm_covariate     = NULL,
                       par_equiclass     = NULL,
                       ...) {

    stopifnot(!is.null(cfg), is.list(cfg))
    cfg_copy <- cfg                    # work on a copy

    # -----------------------------------------------------------------------
    # 1.  Inferences from data_list
    # -----------------------------------------------------------------------
    if (!is.null(data_list)) {
        stopifnot(is.list(data_list))
        tic   <- data_list$time_ind_covariate
        n_tic <- if (is.null(tic)) 0L else as.integer(ncol(tic) - 1L)
        tdc   <- data_list$time_dep_covariate
        n_tdc <- if (is.null(tdc)) 0L else as.integer(ncol(tdc) - 3L)

        cfg_copy$n_time_ind_covariate <- n_tic
        cfg_copy$n_time_dep_covariate <- n_tdc
        cfg_copy$n_covariate          <- n_tic + n_tdc

        cfg_copy$generate_c2p_contact <- if (is.null(data_list$c2p_contact)) 1L else 0L
        cfg_copy$generate_p2p_contact <- if (is.null(data_list$p2p_contact)) 1L else 0L

        if (is.null(data_list$c2p_contact) && is.null(data_list$p2p_contact))
            cfg_copy$common_contact_history_within_community <- 1L
    }

    # -----------------------------------------------------------------------
    # 2.  Generic direct-assignment extras  (...)
    # -----------------------------------------------------------------------
    extra <- list(...)
    for (nm in names(extra))
        cfg_copy[[nm]] <- extra[[nm]]

    # -----------------------------------------------------------------------
    # 3.  Paths
    # -----------------------------------------------------------------------
    if (!is.null(path_in))  cfg_copy$path_in  <- path_in
    if (!is.null(path_out)) cfg_copy$path_out <- path_out

    # -----------------------------------------------------------------------
    # 4.  Incubation period
    # -----------------------------------------------------------------------
    inc_prov <- !c(is.null(min_inc), is.null(max_inc), is.null(prob_inc))
    if (any(inc_prov)) {
        if (!all(inc_prov))
            stop("min_inc, max_inc, and prob_inc must all be provided together.")
        if (!is.list(prob_inc)) prob_inc <- list(prob_inc)
        n <- length(min_inc)
        if (length(max_inc) != n || length(prob_inc) != n)
            stop("min_inc, max_inc, and prob_inc must all have the same length.")
        cfg_copy$n_inc    <- as.integer(n)
        cfg_copy$min_inc  <- as.integer(min_inc)
        cfg_copy$max_inc  <- as.integer(max_inc)
        cfg_copy$prob_inc <- prob_inc
    }

    # -----------------------------------------------------------------------
    # 5.  Infectious period
    # -----------------------------------------------------------------------
    inf_prov <- !c(is.null(lower_inf), is.null(upper_inf), is.null(prob_inf))
    if (any(inf_prov)) {
        if (!all(inf_prov))
            stop("lower_inf, upper_inf, and prob_inf must all be provided together.")
        if (!is.list(prob_inf)) prob_inf <- list(prob_inf)
        n <- length(lower_inf)
        if (length(upper_inf) != n || length(prob_inf) != n)
            stop("lower_inf, upper_inf, and prob_inf must all have the same length.")
        cfg_copy$n_inf     <- as.integer(n)
        cfg_copy$lower_inf <- as.integer(lower_inf)
        cfg_copy$upper_inf <- as.integer(upper_inf)
        cfg_copy$prob_inf  <- prob_inf
    }

    # -----------------------------------------------------------------------
    # 6.  Covariate specifications
    # -----------------------------------------------------------------------

    # Build covariate name → global index lookup from data_list
    .build_lookup <- function(dl) {
        lookup <- integer(0L)
        if (!is.null(dl)) {
            tic_dl  <- dl$time_ind_covariate
            n_tic_l <- if (is.null(tic_dl)) 0L else ncol(tic_dl) - 1L
            if (!is.null(tic_dl) && n_tic_l > 0L)
                lookup <- c(lookup,
                            setNames(seq_len(n_tic_l), names(tic_dl)[-1L]))
            tdc_dl <- dl$time_dep_covariate
            if (!is.null(tdc_dl)) {
                n_tdc_l <- ncol(tdc_dl) - 3L
                if (n_tdc_l > 0L)
                    lookup <- c(lookup,
                                setNames(n_tic_l + seq_len(n_tdc_l),
                                         names(tdc_dl)[-(1:3)]))
            }
        }
        lookup
    }
    cov_lookup <- .build_lookup(data_list)

    .resolve_cov <- function(spec, arg_name) {
        if (is.character(spec)) {
            if (length(cov_lookup) == 0L)
                stop("data_list must be provided to resolve covariate names ",
                     "for '", arg_name, "'.")
            idx <- cov_lookup[spec]
            bad <- spec[is.na(idx)]
            if (length(bad) > 0L)
                stop("Covariate name(s) not found in data_list for '",
                     arg_name, "': ", paste(bad, collapse = ", "))
            as.integer(idx)
        } else {
            as.integer(spec)
        }
    }

    cov_changed <- FALSE

    if (!is.null(c2p_covariate)) {
        cfg_copy$c2p_covariate   <- .resolve_cov(c2p_covariate,   "c2p_covariate")
        cfg_copy$n_c2p_covariate <- length(cfg_copy$c2p_covariate)
        cov_changed <- TRUE
    }
    if (!is.null(sus_p2p_covariate)) {
        cfg_copy$sus_p2p_covariate   <- .resolve_cov(sus_p2p_covariate, "sus_p2p_covariate")
        cfg_copy$n_sus_p2p_covariate <- length(cfg_copy$sus_p2p_covariate)
        cov_changed <- TRUE
    }
    if (!is.null(inf_p2p_covariate)) {
        cfg_copy$inf_p2p_covariate   <- .resolve_cov(inf_p2p_covariate, "inf_p2p_covariate")
        cfg_copy$n_inf_p2p_covariate <- length(cfg_copy$inf_p2p_covariate)
        cov_changed <- TRUE
    }
    if (!is.null(pat_covariate)) {
        cfg_copy$pat_covariate   <- .resolve_cov(pat_covariate, "pat_covariate")
        cfg_copy$n_pat_covariate <- length(cfg_copy$pat_covariate)
        cov_changed <- TRUE
    }
    if (!is.null(imm_covariate)) {
        cfg_copy$imm_covariate   <- .resolve_cov(imm_covariate, "imm_covariate")
        cfg_copy$n_imm_covariate <- length(cfg_copy$imm_covariate)
        cov_changed <- TRUE
    }

    if (cov_changed) {
        cfg_copy$n_p2p_covariate <- as.integer(
            cfg_copy$n_sus_p2p_covariate +
            cfg_copy$n_inf_p2p_covariate +
            cfg_copy$n_int_p2p_covariate)
        cfg_copy$n_par <- as.integer(
            cfg_copy$n_b_mode        + cfg_copy$n_p_mode       +
            cfg_copy$n_u_mode        + cfg_copy$n_q_mode        +
            cfg_copy$n_c2p_covariate + cfg_copy$n_sus_p2p_covariate +
            cfg_copy$n_inf_p2p_covariate + cfg_copy$n_int_p2p_covariate +
            cfg_copy$n_pat_covariate + cfg_copy$n_imm_covariate)
    }

    # -----------------------------------------------------------------------
    # 7.  par_equiclass  — members may be integer indices or parameter labels
    # -----------------------------------------------------------------------
    if (!is.null(par_equiclass)) {

        # Build full parameter label vector (one entry per sequential index)
        .lbl <- function(prefix, n) {
            if (n == 0L) character(0L)
            else if (n == 1L) prefix
            else paste0(prefix, seq_len(n))
        }

        # Helper: global covariate index → label (uses data_list if available)
        n_tic_cfg <- cfg_copy$n_time_ind_covariate
        .cov_lbl <- function(k) {
            if (!is.null(data_list)) {
                if (k <= n_tic_cfg) {
                    tic <- data_list$time_ind_covariate
                    if (!is.null(tic) && ncol(tic) >= k + 1L)
                        return(names(tic)[k + 1L])
                } else {
                    tdc <- data_list$time_dep_covariate
                    col <- k - n_tic_cfg + 3L
                    if (!is.null(tdc) && ncol(tdc) >= col)
                        return(names(tdc)[col])
                }
            }
            paste0("x", k)
        }
        .vcov <- function(ids, n) {
            if (n == 0L) character(0L)
            else vapply(ids[seq_len(n)], .cov_lbl, character(1L))
        }

        full_par_labels <- c(
            .lbl("b",   cfg_copy$n_b_mode),
            .lbl("p",   cfg_copy$n_p_mode),
            .lbl("u",   cfg_copy$n_u_mode),
            .lbl("q",   cfg_copy$n_q_mode),
            if (cfg_copy$n_c2p_covariate     > 0L)
                paste0("c2p_",
                       .vcov(cfg_copy$c2p_covariate, cfg_copy$n_c2p_covariate)),
            if (cfg_copy$n_sus_p2p_covariate > 0L)
                paste0("p2p_s_",
                       .vcov(cfg_copy$sus_p2p_covariate, cfg_copy$n_sus_p2p_covariate)),
            if (cfg_copy$n_inf_p2p_covariate > 0L)
                paste0("p2p_i_",
                       .vcov(cfg_copy$inf_p2p_covariate, cfg_copy$n_inf_p2p_covariate)),
            if (!is.null(cfg_copy$interaction) && cfg_copy$n_int_p2p_covariate > 0L)
                vapply(seq_len(cfg_copy$n_int_p2p_covariate), function(j)
                    paste0("p2p_int_",
                           .cov_lbl(cfg_copy$interaction[[j]][1L]),
                           "_",
                           .cov_lbl(cfg_copy$interaction[[j]][2L])), character(1L)),
            .lbl("pat", cfg_copy$n_pat_covariate),
            .lbl("imm", cfg_copy$n_imm_covariate)
        )

        lbl2idx <- setNames(seq_along(full_par_labels), full_par_labels)

        resolved <- lapply(par_equiclass, function(ec) {
            m <- ec$member
            if (is.character(m)) {
                idx <- lbl2idx[m]
                bad <- m[is.na(idx)]
                if (length(bad) > 0L)
                    stop("Parameter label(s) not found for par_equiclass: ",
                         paste(bad, collapse = ", "),
                         "\nAvailable labels: ",
                         paste(full_par_labels, collapse = ", "))
                m <- as.integer(idx)
            } else {
                m <- as.integer(m)
            }
            list(size = length(m), member = m)
        })

        cfg_copy$par_equiclass   <- resolved
        cfg_copy$n_par_equiclass <- as.integer(length(resolved))
    }

    # -----------------------------------------------------------------------
    # 8.  Refresh par_labels when data_list is available
    # -----------------------------------------------------------------------
    if (!is.null(data_list))
        cfg_copy <- update_par_labels(cfg_copy, data_list)

    cfg_copy
}
