#' Generate a Pseudo-population for Chain Binomial Simulation
#'
#' Creates a synthetic population of communities, assigns treatment
#' covariates, and seeds index cases. This is the first step in a
#' chain binomial simulation study. The population can then be passed
#' to the TranStat estimation routines.
#'
#' @param n_community Integer. Number of independent communities (e.g., households).
#' @param community_size Integer. Number of individuals per community.
#' @param day_epi_stop Integer. Last day of the epidemic follow-up period.
#' @param case_ascertained Integer (0 or 1). If 1, the first member of each
#'   community is pre-assigned as the index case with illness onset on day 1
#'   (case-ascertained design). Default 0.
#' @param cluster_randomization Integer (0 or 1). If 1, all index cases in a
#'   community receive the same treatment assignment, and all contacts receive
#'   a (potentially different) shared assignment. If 0, each individual is
#'   randomized independently. Default 0.
#' @param seed Integer. Random seed. Default 123456789L.
#' @param prop_idx Double. Probability of assigning covariate value 1 to an
#'   index case. Default 0.8.
#' @param prop_contact Double. Probability of assigning covariate value 1 to a
#'   household contact. Default 0.5.
#'
#' @return A named list with four data frames:
#' \describe{
#'   \item{pop}{One row per person. Community assignment and disease-related
#'     information: \code{id}, \code{community}, \code{pre_immune},
#'     \code{infection}, \code{symptom}, \code{day_ill}, \code{exit},
#'     \code{day_exit}, \code{idx}, \code{u_mode}, \code{q_mode},
#'     \code{weight}, \code{ignore}.}
#'   \item{time_ind_covariate}{One row per person. Columns: \code{id},
#'     \code{value}. The covariate is 1 with probability \code{prop_idx}
#'     for index cases and \code{prop_contact} for contacts.}
#'   \item{community}{One row per community. Columns: \code{community},
#'     \code{day_epi_start}, \code{day_epi_stop}, \code{day_last_followup},
#'     \code{c2p_group}.}
#'   \item{time_dep_covariate}{Interval-format time-dependent covariate
#'     (standard normal). Columns: \code{id}, \code{day_start}, \code{day_stop},
#'     \code{value}. Each row records the covariate value over the interval
#'     [\code{day_start}, \code{day_stop}]. In this simulation each interval
#'     spans a single day (\code{day_start == day_stop}), but the format
#'     generalises to covariates that are constant across multi-day segments
#'     (e.g., vaccination status).}
#' }
#'
#' @examples
#' \dontrun{
#' pop <- gen_population(
#'   n_community       = 100,
#'   community_size    = 5,
#'   day_epi_stop      = 14,
#'   case_ascertained  = 1L,
#'   cluster_randomization = 0L
#' )
#' head(pop$pop)
#' }
#'
#' @export
gen_population <- function(n_community, community_size, day_epi_stop,
                           case_ascertained = 0L, cluster_randomization = 0L,
                           seed = 123456789L, prop_idx = 0.8, prop_contact = 0.5) {
    res <- .Call("r_gen_population",
                 as.integer(n_community),
                 as.integer(community_size),
                 as.integer(day_epi_stop),
                 as.integer(case_ascertained),
                 as.integer(cluster_randomization),
                 as.integer(seed),
                 as.double(prop_idx),
                 as.double(prop_contact))

    # Rename the single time-independent covariate column to "x1"
    # and the single time-dependent covariate column to "x2", matching
    # the default sequential naming used by read_population().
    if (!is.null(res$time_ind_covariate))
        names(res$time_ind_covariate) <- c("id", "x1")
    if (!is.null(res$time_dep_covariate))
        names(res$time_dep_covariate) <- c("id", "day_start", "day_stop", "x2")

    res
}
