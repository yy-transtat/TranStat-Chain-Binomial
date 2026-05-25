#' Simulate an Epidemic on a Pseudo-population
#'
#' Writes the pseudo-population produced by \code{\link{gen_population}} to
#' the input directory specified in \code{config_file}, calls the TranStat C
#' simulation engine via \code{\link{transtat}}, then reads back every
#' per-replicate output file written to the output directory
#' (\code{sim_pop_0.txt}, \code{sim_pop_1.txt}, …) and returns them as a
#' single tidy data frame.
#'
#' The number of replicates, epidemic parameters, incubation and infectious
#' period distributions, and covariate structure are all governed by the
#' \code{config_file}.  The file must contain at minimum:
#' \itemize{
#'   \item \code{# perform-simulation} set to \code{1}
#'   \item \code{# output-simulation-data} set to \code{1}
#'   \item \code{# simulation-only} set to \code{1} if estimation should be skipped
#' }
#'
#' @param pop_list Named list returned by \code{\link{gen_population}}.
#' @param config_file Character. Path to the TranStat \code{config.file}.
#' @param serial_number Integer. Run identifier appended to some output file
#'   names when running multiple replicates. Default \code{1L}.
#' @param seed Integer. Random seed passed to the C RNG. Default
#'   \code{12345678L}.
#'
#' @return A data frame with one row per person per simulation replicate.
#'   Columns:
#'   \describe{
#'     \item{\code{sim_id}}{Zero-based replicate index matching the file
#'       suffix in \code{sim_pop_<sim_id>.txt}.}
#'     \item{\code{id}}{Person identifier (0-based).}
#'     \item{\code{community}}{Community (household) identifier.}
#'     \item{\code{idx}}{Index-case flag: 1 if this person is a pre-designated
#'       index case, 0 otherwise.}
#'     \item{\code{infection}}{1 if infected during the simulated epidemic,
#'       0 otherwise.}
#'     \item{\code{day_infection}}{Day of infection (\code{NA} if not
#'       infected).}
#'     \item{\code{symptom}}{1 if the infection is symptomatic, 0 if
#'       asymptomatic or uninfected.}
#'     \item{\code{day_ill}}{Illness / infectiousness onset day (\code{NA}
#'       if not infected).}
#'     \item{\code{x1}, \code{x2}, …}{Time-independent covariate values
#'       at the time of outcome, one column per covariate
#'       (\code{n_time_ind_covariate} columns from the config file).}
#'   }
#'
#' @seealso \code{\link{gen_population}}, \code{\link{write_population}},
#'   \code{\link{read_config}}, \code{\link{transtat}}
#'
#' @examples
#' \dontrun{
#' pop <- gen_population(
#'   n_community      = 100,
#'   community_size   = 5,
#'   day_epi_stop     = 14,
#'   case_ascertained = 1L
#' )
#' cfg_file <- system.file("extdata", "CaseStudy1", "config.file", package = "ChainBinomial")
#' sims <- simulate_epidemic(pop, cfg_file)
#'
#' # Overall infection and secondary attack rate across all replicates
#' contacts <- subset(sims, idx == 0)
#' mean(contacts$infection)
#'
#' # Per-replicate SAR
#' aggregate(infection ~ sim_id, data = contacts, FUN = mean)
#' }
#'
#' @export
simulate_epidemic <- function(pop_list, config_file,
                               serial_number = 1L, seed = 12345678L) {

    cfg <- read_config(config_file)

    # ---- 1.  Write population files to the input directory ----
    write_population(pop_list, dir = cfg$path_in)

    # ---- 2.  Run TranStat (simulation controlled by config_file) ----
    transtat(config_file,
             serial_number = as.integer(serial_number),
             seed          = as.integer(seed))

    # ---- 3.  Locate per-replicate output files ----
    sim_files <- sort(
        list.files(cfg$path_out,
                   pattern  = "^sim_pop_[0-9]+\\.txt$",
                   full.names = TRUE)
    )
    if (length(sim_files) == 0L)
        stop("No sim_pop_*.txt files found in '", cfg$path_out, "'. ",
             "Verify that config.file has 'output-simulation-data = 1' ",
             "and 'perform-simulation = 1'.")

    # ---- 4.  Read results ----
    # Column layout written by core.h:
    #   id  community  idx  infection  day_infection  symptom  day_ill
    #   [time_ind_cov_0  time_ind_cov_1  ...]
    n_tic     <- cfg$n_time_ind_covariate
    col_names <- c("id", "community", "idx",
                   "infection", "day_infection",
                   "symptom",   "day_ill",
                   if (n_tic > 0L) paste0("x", seq_len(n_tic)) else character(0))

    sims <- lapply(seq_along(sim_files), function(k) {
        df        <- read.table(sim_files[k], header = FALSE,
                                col.names = col_names)
        df$sim_id <- k - 1L

        # Replace MISSING sentinel (99999) with NA
        df$day_infection[df$day_infection == 99999L] <- NA_integer_
        df$day_ill[df$day_ill       == 99999L] <- NA_integer_
        df
    })

    out          <- do.call(rbind, sims)
    rownames(out) <- NULL
    out[, c("sim_id", col_names)]
}
