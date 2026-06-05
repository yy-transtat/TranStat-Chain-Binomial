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
#' @param int_p2p_covariate List of interaction pairs.  Each element is a
#'   length-2 vector \code{c(sus_idx, inf_idx)} giving the global covariate
#'   indices (the sequential numbering across all time-independent then
#'   time-dependent covariates) of the susceptibility and infectiousness
#'   covariates forming the interaction.  Integer form: values are used
#'   directly as global indices.  Character form \code{c("sus_name",
#'   "inf_name")}: each name is looked up in \code{covariate_labels} and
#'   resolved to its global index.
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
                       int_p2p_covariate = NULL,
                       par_equiclass     = NULL,
                       ...) {

    stopifnot(!is.null(cfg), is.list(cfg))
    cfg_copy <- cfg                    # work on a copy

    # Remove *_names fields — these are not stored in any C data structure and
    # are superseded by the auto-generated "class1", "class2", ... labels that
    # write_config() emits.
    cfg_copy$sim_par_names           <- NULL
    cfg_copy$converge_criteria_names <- NULL
    cfg_copy$ini_par_names           <- NULL
    cfg_copy$search_bound_names      <- NULL

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

    # Build covariate name → global index lookup.
    # Priority: (1) data_list column names, (2) cfg_copy$covariate_labels.
    .build_lookup <- function(dl, cfg_arg) {
        if (!is.null(dl)) {
            # Build from data_list column names
            lookup  <- integer(0L)
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
            return(lookup)
        }
        # Fall back to cfg$covariate_labels when data_list is absent
        lbl <- cfg_arg$covariate_labels
        if (!is.null(lbl) && length(lbl) > 0L)
            return(setNames(seq_along(lbl), lbl))
        integer(0L)
    }
    cov_lookup <- .build_lookup(data_list, cfg_copy)

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

    # Interaction pairs: each element is c(sus_idx, inf_idx) where both values
    # are global covariate indices (the same sequential numbering used for
    # time-independent and time-dependent covariates across the whole dataset).
    # Character form: each name is looked up in cov_lookup (name → global index).
    # Integer form: values are already global indices — copy directly.
    if (!is.null(int_p2p_covariate)) {
        if (!is.list(int_p2p_covariate))
            stop("int_p2p_covariate must be a list of length-2 vectors.")
        resolved_int <- lapply(seq_along(int_p2p_covariate), function(k) {
            pair <- int_p2p_covariate[[k]]
            if (length(pair) != 2L)
                stop("int_p2p_covariate element ", k, " must have length 2.")
            if (is.character(pair)) {
                idx <- cov_lookup[pair]
                bad <- pair[is.na(idx)]
                if (length(bad) > 0L)
                    stop("Covariate name(s) not found for int_p2p_covariate ",
                         "element ", k, ": ", paste(bad, collapse = ", "))
                as.integer(idx)
            } else {
                as.integer(pair)   # already global covariate indices
            }
        })
        cfg_copy$interaction         <- resolved_int
        cfg_copy$n_int_p2p_covariate <- as.integer(length(resolved_int))
        cov_changed <- TRUE
    }

    if (cov_changed) {
        cfg_copy$n_p2p_covariate <- as.integer(
            cfg_copy$n_sus_p2p_covariate +
            cfg_copy$n_inf_p2p_covariate +
            cfg_copy$n_int_p2p_covariate)
        # n_covariate = n_time_ind_covariate + n_time_dep_covariate;
        # it reflects the data layout, not the covariate-effect counts,
        # so it is not recalculated here.
        cfg_copy$n_par <- as.integer(
            cfg_copy$n_b_mode        + cfg_copy$n_p_mode       +
            cfg_copy$n_u_mode        + cfg_copy$n_q_mode        +
            cfg_copy$n_c2p_covariate + cfg_copy$n_sus_p2p_covariate +
            cfg_copy$n_inf_p2p_covariate + cfg_copy$n_int_p2p_covariate +
            cfg_copy$n_pat_covariate + cfg_copy$n_imm_covariate)
    }

    # When any covariate specification changes the equivalence classes are
    # invalidated — the caller must supply a new par_equiclass.
    if (cov_changed && is.null(par_equiclass))
        stop("par_equiclass must be specified when any covariate specification is ",
             "changed (c2p_covariate, sus_p2p_covariate, inf_p2p_covariate, ",
             "pat_covariate, imm_covariate, or int_p2p_covariate).")

    # -----------------------------------------------------------------------
    # 7.  par_equiclass  — members may be integer indices or parameter labels
    # -----------------------------------------------------------------------

    # Refresh covariate_labels and par_labels now, so that any covariate-spec
    # changes made in sections 1-6 are reflected before we try to resolve
    # character parameter labels in par_equiclass.
    # update_var_par_labels() is safe to call here regardless of whether
    # data_list is NULL: it preserves existing covariate_labels when
    # data_list = NULL and they already have the correct length.
    cfg_copy <- update_var_par_labels(cfg_copy, data_list)

    if (!is.null(par_equiclass)) {

        # cfg_copy$par_labels is current — use it directly as the lookup table.
        full_par_labels <- cfg_copy$par_labels
        lbl2idx         <- setNames(seq_along(full_par_labels), full_par_labels)

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
    # 8.  Final label refresh (safety net — idempotent given the call above).
    # -----------------------------------------------------------------------
    cfg_copy <- update_var_par_labels(cfg_copy, data_list)

    # -----------------------------------------------------------------------
    # 9.  Validate dependent fields against n_par_equiclass
    # -----------------------------------------------------------------------
    n_ec <- cfg_copy$n_par_equiclass
    if (!is.null(n_ec) && n_ec > 0L) {

        is_prob <- .is_prob_equiclass(cfg_copy)

        # sim_par_effective: required when simulation == 1
        if (isTRUE(cfg_copy$simulation == 1L)) {
            sv <- cfg_copy$sim_par_effective
            if (is.null(sv) || length(sv) != n_ec)
                stop("sim_par_effective must be a numeric vector of length ",
                     n_ec, " (n_par_equiclass) when simulation == 1.")
            .check_par_values(as.double(sv), is_prob, "sim_par_effective")
        }

        # ini_par_effective: normalize list → matrix when ini_par_provided == 1
        if (isTRUE(cfg_copy$ini_par_provided == 1L)) {
            ip <- cfg_copy$ini_par_effective
            if (is.null(ip))
                stop("ini_par_effective must be provided when ini_par_provided == 1.")
            if (is.list(ip)) {
                for (k in seq_along(ip))
                    if (length(ip[[k]]) != n_ec)
                        stop("ini_par_effective[[", k, "]] must have length ",
                             n_ec, " (n_par_equiclass).")
                ip <- do.call(rbind, lapply(ip, as.double))
                cfg_copy$ini_par_effective <- ip
                cfg_copy$n_ini             <- as.integer(nrow(ip))
            } else if (is.numeric(ip)) {
                if (is.vector(ip)) {
                    if (length(ip) != n_ec)
                        stop("ini_par_effective must have length ",
                             n_ec, " (n_par_equiclass).")
                    ip <- matrix(as.double(ip), nrow = 1L)
                    cfg_copy$ini_par_effective <- ip
                    cfg_copy$n_ini             <- 1L
                } else if (is.matrix(ip)) {
                    if (ncol(ip) != n_ec)
                        stop("ini_par_effective must have ", n_ec,
                             " columns (n_par_equiclass).")
                }
            } else {
                stop("ini_par_effective must be a numeric vector, matrix, or ",
                     "list of numeric vectors.")
            }
            for (i in seq_len(nrow(cfg_copy$ini_par_effective)))
                .check_par_values(cfg_copy$ini_par_effective[i, ], is_prob,
                                  paste0("ini_par_effective row ", i))
        } else {
            # ini_par_provided == 0: set to NULL
            cfg_copy$ini_par_effective <- NULL
            cfg_copy$n_ini             <- 1L
        }

        # converge_criteria: required when converge_criteria_provided == 1
        if (isTRUE(cfg_copy$converge_criteria_provided == 1L)) {
            cc <- cfg_copy$converge_criteria
            if (is.null(cc) || length(cc) != n_ec)
                stop("converge_criteria must be a numeric vector of length ",
                     n_ec, " when converge_criteria_provided == 1.")
        }

        # lower/upper_search_bound: required when search_bound_provided == 1
        if (isTRUE(cfg_copy$search_bound_provided == 1L)) {
            lb <- cfg_copy$lower_search_bound
            ub <- cfg_copy$upper_search_bound
            if (is.null(lb) || length(lb) != n_ec)
                stop("lower_search_bound must be a numeric vector of length ",
                     n_ec, " when search_bound_provided == 1.")
            if (is.null(ub) || length(ub) != n_ec)
                stop("upper_search_bound must be a numeric vector of length ",
                     n_ec, " when search_bound_provided == 1.")
            .check_par_values(as.double(lb), is_prob, "lower_search_bound")
            .check_par_values(as.double(ub), is_prob, "upper_search_bound")
        }
    }

    cfg_copy
}

# print covariate or parameter labels for given indices
pr.label <- function(id, labels){
  ifelse(length(id)==1, labels[id],
         paste(labels[id], collapse = ', '))
}
# a function displaying defined types of covariates and the numbers of these covariates
show_cfg_covariates <- function(cfg) {
  data.frame(Type_Covariates=
               c('time-independent', 'time-dependent',
                 'c2p', 'p2p susceptibility',
                 'p2p infectivity', 'p2p', 'pathogenicity',
                 'preseason immunity'),
             Size=
               c(cfg$n_time_ind_covariate, cfg$n_time_dep_covariate,
                 cfg$n_c2p_covariate, cfg$n_sus_p2p_covariate,
                 cfg$n_inf_p2p_covariate, cfg$n_p2p_covariate,
                 cfg$n_pat_covariate, cfg$n_imm_covariate),
             Variables=
               c('', '',
                 paste(cfg$c2p_covariate, collapse = ','),
                 paste(cfg$sus_p2p_covariate, collapse = ','),
                 paste(cfg$inf_p2p_covariate, collapse = ','),
                 '',
                 paste(cfg$pat_covariate, collapse = ','),
                 paste(cfg$imm_covariate, collapse = ',')),
             Lables=
               c('', '',
                 pr.label(cfg$c2p_covariate, cfg$covariate_labels),
                 pr.label(cfg$sus_p2p_covariate, cfg$covariate_labels),
                 pr.label(cfg$inf_p2p_covariate, cfg$covariate_labels),
                 '',
                 pr.label(cfg$pat_covariate, cfg$covariate_labels),
                 pr.label(cfg$imm_covariate, cfg$covariate_labels)))
}

# a function displaying equivalence classes in a more friendly way
show_cfg_par_equiclass <- function(cfg) {
  size <- sapply(cfg$par_equiclass, "[[", "size")
  mem.lst <- sapply(cfg$par_equiclass, "[[", "member")
  members <- sapply(mem.lst, paste, collapse = ", ")
  label <- sapply(mem.lst, pr.label, labels=cfg$par_labels)
  data.frame(size, members, label)
}

# ---------------------------------------------------------------------------
# Internal helpers for parameter-type classification and value validation.
# Used by both set_config() and read_config() — not exported.
# ---------------------------------------------------------------------------

# Returns a logical vector of length n_par_equiclass:
#   TRUE  = probability parameter (label matches "^[bpuq][0-9]*$") → (0, 1)
#   FALSE = odds-ratio / covariate parameter → > 0
.is_prob_equiclass <- function(cfg) {
    n_ec <- if (is.null(cfg$n_par_equiclass)) 0L else cfg$n_par_equiclass
    if (is.null(cfg$par_equiclass) || n_ec == 0L)
        return(logical(0L))
    par_lbl <- cfg$par_labels
    vapply(cfg$par_equiclass, function(ec) {
        lbl <- if (!is.null(par_lbl) && length(par_lbl) >= ec$member[1L])
                   par_lbl[ec$member[1L]]
               else ""
        grepl("^[bpuq][0-9]*$", lbl)
    }, logical(1L))
}

# Check that each value in vals satisfies the constraint implied by is_prob.
# vals    : numeric vector, one value per equivalence class
# is_prob : logical vector (TRUE = probability in (0,1); FALSE = OR > 0)
# what    : label for error messages (e.g. "ini_par_effective row 1")
.check_par_values <- function(vals, is_prob, what) {
    for (k in seq_along(vals)) {
        v <- vals[k]
        if (is.na(v)) next
        if (is_prob[k]) {
            if (v <= 0 || v >= 1)
                stop(what, " class ", k, ": value ", v,
                     " is not in (0, 1) — expected a probability.")
        } else {
            if (v <= 0)
                stop(what, " class ", k, ": value ", v,
                     " is not > 0 — expected an odds ratio.")
        }
    }
    invisible(NULL)
}
