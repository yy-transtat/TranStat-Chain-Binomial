#' Update Covariate and Parameter Labels Using Data Column Names
#'
#' Rebuilds two R-only fields in \code{cfg}:
#' \describe{
#'   \item{\code{cfg$covariate_labels}}{Character vector of length
#'     \code{n_covariate}.  Entry \eqn{k} is the name for global covariate
#'     index \eqn{k}: time-independent covariates occupy positions
#'     \eqn{1, \ldots, n_{tic}} (in the order they appear as columns of
#'     \code{data_list$time_ind_covariate}), followed by time-dependent
#'     covariates at positions \eqn{n_{tic}+1, \ldots, n_{tic}+n_{tdc}}.
#'     When \code{data_list} is \code{NULL} or a column is unavailable the
#'     placeholder \code{"x\{k\}"} is used.}
#'   \item{\code{cfg$par_labels}}{Character vector of length \code{n_par}.
#'     Labels for every model parameter in canonical order (b / p / u / q /
#'     c2p / p2p_s / p2p_i / p2p_int / pat / imm), derived from
#'     \code{covariate_labels}.  Base parameters use the option-2 convention:
#'     the trailing digit is omitted when there is only one mode (e.g.\ "b"
#'     not "b1").}
#' }
#'
#' \code{\link{read_config}} calls this function internally with
#' \code{data_list = NULL} to set placeholder labels.  Once
#' \code{\link{read_population}} has been called and the data frames carry
#' meaningful column names, call this function again to replace the
#' placeholders:
#'
#' \preformatted{
#' cfg       <- read_config(config_file)
#' data_list <- read_population(cfg, names_tid = c("age", "sex"),
#'                                   names_tdp = c("antiviral"))
#' cfg       <- update_var_par_labels(cfg, data_list)
#' cfg$covariate_labels   # "age"  "sex"  "antiviral"
#' cfg$par_labels         # "b"  "p"  "c2p_age"  "p2p_s_sex"  "p2p_i_antiviral"
#' }
#'
#' \code{\link{ChainBinomial}} calls this function automatically.  Users who
#' work directly with \code{\link{estimate_single}} or
#' \code{\link{simulate_epidemics}} should call it manually after
#' \code{\link{read_population}}.
#'
#' @param cfg Named list returned by \code{\link{read_config}}.
#' @param data_list Named list returned by \code{\link{read_population}} or
#'   \code{\link{gen_population}}, or \code{NULL} (default).  When
#'   \code{NULL}, placeholder labels \code{"x1"}, \code{"x2"}, \ldots are
#'   used for all covariates.
#'
#' @return The \code{cfg} list with \code{cfg$covariate_labels} and
#'   \code{cfg$par_labels} updated.  Assign the return value:
#'   \code{cfg <- update_var_par_labels(cfg, data_list)}.
#'
#' @seealso \code{\link{read_config}}, \code{\link{read_population}},
#'   \code{\link{ChainBinomial}}, \code{\link{estimate_single}},
#'   \code{\link{set_config}}
#'
#' @export
update_var_par_labels <- function(cfg, data_list = NULL) {

    stopifnot(!is.null(cfg), is.list(cfg))
    if (!is.null(data_list)) stopifnot(is.list(data_list))

    # Helper: option-2 base-parameter label vector.
    .lbl <- function(prefix, n) {
        if (n == 0L) character(0L)
        else if (n == 1L) prefix
        else paste0(prefix, seq_len(n))
    }

    n_tic <- cfg$n_time_ind_covariate
    n_tdc <- cfg$n_time_dep_covariate
    n_cov <- n_tic + n_tdc

    # -------------------------------------------------------------------
    # 1.  Build covariate_labels (one entry per global covariate index)
    # -------------------------------------------------------------------
    cov_lbl <- character(n_cov)

    tic <- if (!is.null(data_list)) data_list$time_ind_covariate else NULL
    for (k in seq_len(n_tic)) {
        cov_lbl[k] <-
            if (!is.null(tic) && ncol(tic) >= k + 1L) names(tic)[k + 1L]
            else paste0("x", k)
    }

    tdc <- if (!is.null(data_list)) data_list$time_dep_covariate else NULL
    for (k in seq_len(n_tdc)) {
        col <- k + 3L   # columns 1-3 are id / day_start / day_stop
        cov_lbl[n_tic + k] <-
            if (!is.null(tdc) && ncol(tdc) >= col) names(tdc)[col]
            else paste0("x", n_tic + k)
    }

    cfg$covariate_labels <- cov_lbl

    # -------------------------------------------------------------------
    # 2.  Build par_labels using covariate_labels for lookup
    # -------------------------------------------------------------------
    cov_label <- function(k) cov_lbl[k]   # global index → label

    cfg$par_labels <- c(
        .lbl("b", cfg$n_b_mode),
        .lbl("p", cfg$n_p_mode),
        .lbl("u", cfg$n_u_mode),
        .lbl("q", cfg$n_q_mode),
        if (cfg$n_c2p_covariate     > 0L)
            paste0("c2p_",
                   vapply(cfg$c2p_covariate[seq_len(cfg$n_c2p_covariate)],
                          cov_label, character(1L))),
        if (cfg$n_sus_p2p_covariate > 0L)
            paste0("p2p_s_",
                   vapply(cfg$sus_p2p_covariate[seq_len(cfg$n_sus_p2p_covariate)],
                          cov_label, character(1L))),
        if (cfg$n_inf_p2p_covariate > 0L)
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

    cfg
}
