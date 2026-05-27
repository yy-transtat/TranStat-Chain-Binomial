#' Simulate Epidemics and Estimate Parameters
#'
#' Runs \code{n_sim} independent chain-binomial epidemics on the supplied
#' population, estimates transmission parameters from each simulated dataset,
#' and returns the results in tidy wide-format data frames suitable for
#' simulation studies and empirical power calculations.
#'
#' For each iteration \code{k} the function:
#' \enumerate{
#'   \item Calls \code{\link{simulate_single}} with seed \code{seed + k - 1}.
#'   \item Passes the simulated outcomes (together with any contact histories
#'     carried in \code{data_list}) to \code{\link{estimate_single}} using the
#'     same seed.
#' }
#'
#' @param data_list Named list as returned by \code{\link{gen_population}} or
#'   \code{\link{read_population}}.  Elements \code{$c2p_contact} and
#'   \code{$p2p_contact}, when non-\code{NULL}, are forwarded to both
#'   \code{simulate_single} and \code{estimate_single} so that the same
#'   contact structure is used for simulation and estimation.
#' @param cfg Named list returned by \code{\link{read_config}}.
#' @param i_inc Integer (1-based). Incubation-period group index. Default
#'   \code{1L}.
#' @param i_inf Integer (1-based). Infectious-period group index. Default
#'   \code{1L}.
#' @param seed Integer. Base RNG seed. Iteration \code{k} uses
#'   \code{seed + k - 1}. Default \code{12345678L}.
#' @param n_sim Integer. Number of simulations. When \code{NULL} (the
#'   default) the value is taken from \code{cfg$n_simulation}.
#'
#' @return A named list with three data frames, each with one row per
#'   simulation and \code{i_inc}, \code{i_inf}, \code{iter} as the first
#'   three columns.
#' \describe{
#'   \item{\code{estimates}}{Wide data frame of parameter estimates.  After
#'     the index columns the parameters appear in order, each contributing
#'     four consecutive columns: \code{<par>.estimate}, \code{<par>.se},
#'     \code{<par>.ci_lower}, \code{<par>.ci_upper}.}
#'   \item{\code{var}}{Covariance matrix on the probability / odds-ratio
#'     scale, stored row-by-row.  After the index columns the entries are
#'     named \code{<par_i>.<par_j>}.}
#'   \item{\code{var_logit}}{Same layout as \code{var} but on the logit /
#'     log scale.}
#' }
#'
#' @seealso \code{\link{simulate_single}}, \code{\link{estimate_single}},
#'   \code{\link{gen_population}}, \code{\link{read_config}}
#'
#' @examples
#' \dontrun{
#' cfg <- read_config(system.file("extdata", "CaseStudy1", "config.file",
#'                                package = "ChainBinomial"))
#' dat <- gen_population(n_community = 100, community_size = 5,
#'                       day_epi_stop = 14, case_ascertained = 1L)
#' out <- simulate_epidemics(dat, cfg, n_sim = 50L, seed = 1L)
#'
#' # Distribution of p1 estimates across simulations
#' hist(out$estimates$p1.estimate)
#'
#' # Mean and SD of p1
#' mean(out$estimates$p1.estimate)
#' sd(out$estimates$p1.estimate)
#' }
#'
#' @export
simulate_epidemics <- function(data_list, cfg,
                                i_inc = 1L, i_inf = 1L,
                                seed  = 12345678L,
                                n_sim = NULL) {

    stopifnot(is.list(data_list), is.list(cfg))
    stopifnot(i_inc >= 1L, i_inc <= cfg$n_inc)
    stopifnot(i_inf >= 1L, i_inf <= cfg$n_inf)

    if (is.null(n_sim)) n_sim <- cfg$n_simulation
    n_sim <- as.integer(n_sim)
    stopifnot(!is.na(n_sim), n_sim >= 1L)

    # Helper: flatten a square matrix row-by-row into a named one-row data frame
    # prefixed with i_inc, i_inf, iter.
    flatten_mat <- function(mat, ii, jj, k) {
        pn   <- rownames(mat)
        nms  <- as.vector(outer(pn, pn, paste, sep = "."))
        vals <- as.vector(t(mat))
        df   <- as.data.frame(matrix(vals, nrow = 1L,
                                     dimnames = list(NULL, nms)))
        cbind(data.frame(i_inc = ii, i_inf = jj, iter = k,
                         stringsAsFactors = FALSE), df)
    }

    res_estimates <- vector("list", n_sim)
    res_var       <- vector("list", n_sim)
    res_var_logit <- vector("list", n_sim)

    for (k in seq_len(n_sim)) {
        iter_seed <- seed + (k - 1L)

        # --- 1. Simulate one epidemic ---
        sim <- simulate_single(data_list, cfg,
                               i_inc = i_inc, i_inf = i_inf,
                               seed  = iter_seed)

        # --- 2. Estimate from simulated outcomes ---
        # Carry contact histories from data_list into the simulated data list
        # so that estimation uses the same contact structure as simulation.
        sim$c2p_contact <- data_list$c2p_contact
        sim$p2p_contact <- data_list$p2p_contact

        fit <- estimate_single(sim, cfg,
                               i_inc = i_inc, i_inf = i_inf,
                               seed  = iter_seed)

        # --- 3. Wide estimates row ---
        # Columns: i_inc, i_inf, iter,
        #          par1.estimate, par1.se, par1.ci_lower, par1.ci_upper,
        #          par2.estimate, ...
        est_df    <- fit$estimates
        par_names <- est_df$parameter
        wide_vals <- unlist(lapply(seq_len(nrow(est_df)), function(j) {
            p <- par_names[j]
            c(setNames(est_df$estimate[j],  paste0(p, ".estimate")),
              setNames(est_df$se[j],         paste0(p, ".se")),
              setNames(est_df$ci_lower[j],   paste0(p, ".ci_lower")),
              setNames(est_df$ci_upper[j],   paste0(p, ".ci_upper")))
        }))
        res_estimates[[k]] <- cbind(
            data.frame(i_inc = i_inc, i_inf = i_inf, iter = k,
                       stringsAsFactors = FALSE),
            as.data.frame(t(wide_vals), stringsAsFactors = FALSE)
        )

        # --- 4. Variance matrices ---
        res_var[[k]]       <- flatten_mat(fit$var,       i_inc, i_inf, k)
        res_var_logit[[k]] <- flatten_mat(fit$var_logit, i_inc, i_inf, k)
    }

    bind_df <- function(lst) {
        result <- do.call(rbind, lst)
        row.names(result) <- NULL
        result
    }

    list(
        estimates = bind_df(res_estimates),
        var       = bind_df(res_var),
        var_logit = bind_df(res_var_logit)
    )
}
