#' Chain-Binomial Transmission Parameter Estimation
#'
#' Runs \code{\link{estimate_single}} for every combination of incubation-period
#' and infectious-period settings declared in the config, and stacks the results
#' into tidy data frames with two additional index columns, \code{i_inc} and
#' \code{i_inf}, identifying the setting used.
#'
#' When the config declares \code{n_inc} incubation-period groups and
#' \code{n_inf} infectious-period groups, \code{n_inc × n_inf} estimation
#' runs are performed.  The inner loop varies \code{i_inf}; the outer loop
#' varies \code{i_inc}.
#'
#' A typical workflow:
#' \preformatted{
#' cfg       <- read_config(config_file)
#' data_list <- read_population(cfg, names_tid = ..., names_tdp = ...)
#' cfg       <- update_var_par_labels(cfg, data_list)
#' out       <- ChainBinomial(data_list, cfg)
#' }
#'
#' @param data_list Named list returned by \code{\link{read_population}} or
#'   \code{\link{gen_population}}.
#' @param cfg Named list returned by \code{\link{read_config}}, optionally
#'   updated with \code{\link{update_var_par_labels}}.
#' @param seed Integer. Random seed for Monte Carlo / MCEM sampling.
#'   Default \code{12345678L}.
#'
#' @return A named list with the following elements.  Every data frame has
#'   \code{i_inc} and \code{i_inf} as its first two columns.
#' \describe{
#'   \item{\code{estimates}}{Data frame combining parameter estimates across
#'     all (i_inc, i_inf) combinations (columns: \code{i_inc}, \code{i_inf},
#'     \code{parameter}, \code{estimate}, \code{se}, \code{ci_lower},
#'     \code{ci_upper}, \code{z}, \code{p_value}).}
#'   \item{\code{SAR}}{Unadjusted secondary attack rates (\code{i_inc},
#'     \code{i_inf}, \code{group}, \code{SAR}, \code{se}, \code{ci_lower},
#'     \code{ci_upper}), or \code{NULL} when no p2p transmission parameters
#'     are present.}
#'   \item{\code{SAR_adjusted}}{Covariate-adjusted SARs (\code{i_inc},
#'     \code{i_inf}, \code{covariate_set}, \code{group}, \code{SAR},
#'     \code{se}, \code{ci_lower}, \code{ci_upper}), or \code{NULL}.}
#'   \item{\code{R0}}{Unadjusted basic reproduction numbers (\code{i_inc},
#'     \code{i_inf}, \code{R0}, \code{se}, \code{ci_lower}, \code{ci_upper}),
#'     or \code{NULL} when no R0 multiplier is provided.}
#'   \item{\code{R0_adjusted}}{Covariate-adjusted R0 values (\code{i_inc},
#'     \code{i_inf}, \code{covariate_set}, \code{R0}, \code{se},
#'     \code{ci_lower}, \code{ci_upper}), or \code{NULL}.}
#'   \item{\code{log_likelihood}}{Data frame with one row per (i_inc, i_inf)
#'     combination (\code{i_inc}, \code{i_inf}, \code{log_likelihood},
#'     \code{error_code}).}
#'   \item{\code{var}}{Data frame with one row per combination. After
#'     \code{i_inc} and \code{i_inf}, the remaining columns are the
#'     covariance matrix on the probability / odds-ratio scale, stored
#'     row-by-row and named \code{par_i.par_j}.}
#'   \item{\code{var_logit}}{Same layout as \code{var} but for the
#'     covariance matrix on the logit / log scale.}
#' }
#'
#' @seealso \code{\link{estimate_single}}, \code{\link{read_config}},
#'   \code{\link{read_population}}, \code{\link{update_var_par_labels}}
#'
#' @examples
#' \dontrun{
#' cfg_file  <- system.file("extdata", "CaseStudy2", "config.file",
#'                          package = "ChainBinomial")
#' cfg       <- read_config(cfg_file)
#' data_list <- read_population(cfg)
#' cfg       <- update_var_par_labels(cfg, data_list)
#' out       <- ChainBinomial(data_list, cfg)
#' out$estimates
#' out$SAR
#' out$R0
#' out$log_likelihood
#' }
#'
#' @export
ChainBinomial <- function(data_list, cfg, seed = 12345678L) {

    stopifnot(is.list(data_list), is.list(cfg))

    n_inc <- cfg$n_inc
    n_inf <- cfg$n_inf
    n_run <- n_inc * n_inf

    # Pre-allocate result containers (one slot per run)
    res_estimates    <- vector("list", n_run)
    res_SAR          <- vector("list", n_run)
    res_SAR_adjusted <- vector("list", n_run)
    res_R0           <- vector("list", n_run)
    res_R0_adjusted  <- vector("list", n_run)
    res_ll           <- vector("list", n_run)
    res_var          <- vector("list", n_run)
    res_var_logit    <- vector("list", n_run)

    # Helper: prepend i_inc / i_inf columns and reset row names
    prepend_idx <- function(df, ii, jj) {
        if (is.null(df)) return(NULL)
        row.names(df) <- NULL
        cbind(i_inc = ii, i_inf = jj, df, stringsAsFactors = FALSE)
    }

    # Helper: flatten a square matrix row-by-row into a one-row data frame.
    # Column names are "row_par.col_par" (e.g. "p1.OR1").
    flatten_mat <- function(mat, ii, jj) {
        pn   <- rownames(mat)                          # parameter names
        nms  <- as.vector(outer(pn, pn, paste, sep = "."))   # row-major grid
        vals <- as.vector(t(mat))                      # row-major values
        df   <- as.data.frame(matrix(vals, nrow = 1L,
                                     dimnames = list(NULL, nms)))
        cbind(i_inc = ii, i_inf = jj, df, stringsAsFactors = FALSE)
    }

    idx <- 0L
    for (ii in seq_len(n_inc)) {
        for (jj in seq_len(n_inf)) {
            idx <- idx + 1L

            fit <- estimate_single(data_list, cfg,
                                   i_inc = ii,
                                   i_inf = jj,
                                   seed  = seed)

            res_estimates[[idx]]    <- prepend_idx(fit$estimates,    ii, jj)
            res_SAR[[idx]]          <- prepend_idx(fit$SAR,          ii, jj)
            res_SAR_adjusted[[idx]] <- prepend_idx(fit$SAR_adjusted, ii, jj)
            res_R0[[idx]]           <- prepend_idx(fit$R0,           ii, jj)
            res_R0_adjusted[[idx]]  <- prepend_idx(fit$R0_adjusted,  ii, jj)
            res_ll[[idx]] <- data.frame(
                i_inc          = ii,
                i_inf          = jj,
                log_likelihood = fit$log_likelihood,
                error_code     = fit$error_code,
                stringsAsFactors = FALSE
            )
            res_var[[idx]]       <- flatten_mat(fit$var,       ii, jj)
            res_var_logit[[idx]] <- flatten_mat(fit$var_logit, ii, jj)
        }
    }

    # Stack non-NULL results
    bind_results <- function(lst) {
        lst <- Filter(Negate(is.null), lst)
        if (length(lst) == 0L) return(NULL)
        result <- do.call(rbind, lst)
        row.names(result) <- NULL
        result
    }

    list(
        estimates      = bind_results(res_estimates),
        SAR            = bind_results(res_SAR),
        SAR_adjusted   = bind_results(res_SAR_adjusted),
        R0             = bind_results(res_R0),
        R0_adjusted    = bind_results(res_R0_adjusted),
        log_likelihood = bind_results(res_ll),
        var            = bind_results(res_var),
        var_logit      = bind_results(res_var_logit)
    )
}
