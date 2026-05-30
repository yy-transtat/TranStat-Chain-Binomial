#' Update Parameter Labels Using Covariate Names from a Data List
#'
#' Rebuilds \code{cfg$par_labels} so that covariate parameters are labelled
#' with the actual column names carried in \code{data_list} (e.g.\
#' \code{"c2p_age"}, \code{"p2p_i_antiviral"}) instead of the default
#' global-index placeholders set by \code{\link{read_config}} (e.g.\
#' \code{"c2p_x3"}, \code{"p2p_i_x8"}).
#'
#' \code{\link{read_config}} builds \code{cfg$par_labels} before any data
#' file has been read, so covariate entries necessarily carry placeholder names
#' of the form \code{x\{k\}} (global covariate index).  Once
#' \code{\link{read_population}} has been called and the data frames in
#' \code{data_list} carry meaningful column names (set via the
#' \code{names_tid} / \code{names_tdp} arguments), call this function to
#' synchronise \code{cfg$par_labels} with those names.
#'
#' \code{\link{ChainBinomial}} calls this function automatically.  Users who
#' work directly with \code{\link{estimate_single}} or
#' \code{\link{simulate_epidemics}} should call it manually after
#' \code{\link{read_population}}:
#'
#' \preformatted{
#' cfg       <- read_config(config_file)
#' data_list <- read_population(cfg, names_tid = c("age", "sex"),
#'                                   names_tdp = c("antiviral"))
#' cfg       <- update_par_labels(cfg, data_list)
#' cfg$par_labels   # now shows "c2p_age", "p2p_i_antiviral", etc.
#' }
#'
#' @param cfg Named list returned by \code{\link{read_config}}.
#' @param data_list Named list returned by \code{\link{read_population}} or
#'   \code{\link{gen_population}}.  The column names of
#'   \code{data_list$time_ind_covariate} and
#'   \code{data_list$time_dep_covariate} are used to resolve covariate labels.
#'
#' @return The \code{cfg} list with \code{cfg$par_labels} updated in place.
#'   Assign the return value: \code{cfg <- update_par_labels(cfg, data_list)}.
#'
#' @seealso \code{\link{read_config}}, \code{\link{read_population}},
#'   \code{\link{ChainBinomial}}, \code{\link{estimate_single}}
#'
#' @export
update_par_labels <- function(cfg, data_list) {

    stopifnot(is.list(cfg), is.list(data_list))

    # Helper: option-2 base-parameter label vector.
    # Returns the bare prefix for n == 1, numbered labels for n > 1.
    .lbl <- function(prefix, n) {
        if (n == 0L) character(0L)
        else if (n == 1L) prefix
        else paste0(prefix, seq_len(n))
    }

    # Helper: resolve global covariate index k → column name from data_list.
    # Falls back to "x{k}" when the data frame is unavailable.
    .n_tic <- cfg$n_time_ind_covariate
    cov_label <- function(k) {
        if (k <= .n_tic) {
            tic <- data_list$time_ind_covariate
            if (!is.null(tic) && ncol(tic) >= k + 1L)
                return(names(tic)[k + 1L])
        } else {
            tdc <- data_list$time_dep_covariate
            col  <- k - .n_tic + 3L      # cols 1-3 are id/day_start/day_stop
            if (!is.null(tdc) && ncol(tdc) >= col)
                return(names(tdc)[col])
        }
        paste0("x", k)                   # fallback
    }

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
