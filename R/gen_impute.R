#' Generate an Imputation Table for Asymptomatic Infections
#'
#' Builds the \code{impute.dat} data frame required by the TranStat EM
#' algorithm for simulated data that contains asymptomatic infections.  Only
#' asymptomatic secondary cases (\code{infection == 1} and \code{symptom ==
#' 0}) receive an imputation row; index cases and uninfected individuals are
#' excluded.
#'
#' Because the illness-onset day of an asymptomatic infection is not observed,
#' the EM algorithm needs a plausible range of onset days.  This function
#' derives that range from the simulated infection day and the mean of the
#' incubation period distribution for group \code{i_inc}:
#'
#' \deqn{
#'   \textrm{possible\_asym\_start} =
#'     \max(\,\textrm{day\_infection} + \bar{L} - 3,\; \textrm{day\_epi\_start}\,)
#' }
#' \deqn{
#'   \textrm{possible\_asym\_stop} =
#'     \min(\,\textrm{day\_infection} + \bar{L} + 3,\; \textrm{day\_epi\_stop}\,)
#' }
#'
#' where \eqn{\bar{L}} is the mean latent period (rounded to the nearest
#' integer) and \code{day\_epi\_start} / \code{day\_epi\_stop} are the
#' epidemic window of the individual's community.
#'
#' @param pop_list Named list containing at least \code{$pop} and
#'   \code{$community}, as returned by \code{\link{simulate_single}} or
#'   \code{\link{gen_population}}.
#' @param cfg Named list returned by \code{\link{read_config}}.
#' @param i_inc Integer (1-based).  Which incubation-period group to use when
#'   computing the mean latent period.  Default \code{1L}.
#'
#' @return A data frame with nine columns and one row per asymptomatic
#'   infection:
#' \describe{
#'   \item{\code{person_id}}{0-based person identifier.}
#'   \item{\code{possible_pre_immune}}{Always 0.}
#'   \item{\code{possible_escape}}{Always 0.}
#'   \item{\code{possible_sym}}{Always 0.}
#'   \item{\code{possible_sym_start}}{Always 0.}
#'   \item{\code{possible_sym_stop}}{Always 0.}
#'   \item{\code{possible_asym}}{Always 1.}
#'   \item{\code{possible_asym_start}}{Earliest plausible onset day, clamped to
#'     \code{day_epi_start} of the person's community.}
#'   \item{\code{possible_asym_stop}}{Latest plausible onset day, clamped to
#'     \code{day_epi_stop} of the person's community.}
#' }
#' Returns \code{NULL} when there are no asymptomatic infections.
#'
#' @seealso \code{\link{simulate_single}}, \code{\link{write_population}},
#'   \code{\link{read_config}}
#'
#' @examples
#' \dontrun{
#' cfg_file <- system.file("extdata", "CaseStudy1", "config.file",
#'                         package = "ChainBinomial")
#' cfg <- read_config(cfg_file)
#' pop <- gen_population(n_community = 100, community_size = 5,
#'                       day_epi_stop = 14, case_ascertained = 1L)
#' out <- simulate_single(pop, cfg, seed = 42L)
#'
#' imp <- gen_impute(out, cfg, i_inc = 1L)
#' head(imp)
#'
#' # Attach and write — write_population() will write impute.dat automatically
#' out$impute <- imp
#' write_population(out, dir = "data/sim01")
#' }
#'
#' @export
gen_impute <- function(pop_list, cfg, i_inc = 1L) {

    stopifnot(is.list(pop_list), is.list(cfg))
    stopifnot(is.data.frame(pop_list$pop), is.data.frame(pop_list$community))
    stopifnot(i_inc >= 1L, i_inc <= cfg$n_inc)

    pop  <- pop_list$pop
    comm <- pop_list$community

    # ---- Mean latent period for incubation group i_inc ----
    days        <- seq.int(cfg$min_inc[i_inc], cfg$max_inc[i_inc])
    probs       <- cfg$prob_inc[[i_inc]]
    mean_latent <- round(sum(days * probs / sum(probs)))   # nearest integer day

    # ---- Asymptomatic infections only ----
    asym <- pop[pop$infection == 1L & pop$symptom == 0L, ]
    if (nrow(asym) == 0L) return(NULL)

    # ---- Community epidemic window for each asymptomatic person ----
    comm_idx  <- match(asym$community, comm$community)
    epi_start <- comm$day_epi_start[comm_idx]
    epi_stop  <- comm$day_epi_stop[comm_idx]

    # ---- Plausible onset range, clamped to epidemic window ----
    asym_start <- pmax(asym$day_infection + mean_latent - 3L, epi_start+cfg$min_inc[i_inc])
    asym_stop  <- pmin(asym$day_infection + mean_latent + 3L, epi_stop)

    data.frame(
        person_id           = asym$id,
        possible_pre_immune = 0L,
        possible_escape     = 0L,
        possible_sym        = 0L,
        possible_sym_start  = 0L,
        possible_sym_stop   = 0L,
        possible_asym       = 1L,
        possible_asym_start = as.integer(asym_start),
        possible_asym_stop  = as.integer(asym_stop),
        row.names           = NULL
    )
}
