#' Read a TranStat Population Data Set into R
#'
#' Parses all seven TranStat input files found in the \code{path_in} directory
#' of \code{cfg} and returns them as a named list of data frames.  The
#' function mirrors \code{core.h} lines 201–1083 of the TranStat C engine but
#' stops after reading — no C structures are built and no validation beyond
#' file existence is performed.
#'
#' Missing files are silently returned as \code{NULL}.  This is the intended
#' behaviour for optional files (\code{c2p_contact.dat},
#' \code{p2p_contact.dat}, \code{impute.dat}) when the corresponding contact
#' history or imputation options are auto-generated or unused.
#'
#' @section File formats:
#' All files are whitespace-delimited with no header line.
#'
#' \strong{pop.dat} — one row per individual:
#' \code{id community pre_immune infection symptom day_ill exit day_exit idx
#' u_mode q_mode weight ignore}
#'
#' \strong{community.dat} — one row per community:
#' \code{community day_epi_start day_epi_stop day_last_followup c2p_group}
#'
#' \strong{time_ind_covariate.dat} — one row per individual:
#' \code{id value} (single covariate) or
#' \code{id value_1 value_2 …} (multiple covariates)
#'
#' \strong{time_dep_covariate.dat} — one row per individual per time segment:
#' \code{id day_start day_stop value} (single covariate) or
#' \code{id day_start day_stop value_1 value_2 …} (multiple covariates)
#'
#' \strong{c2p_contact.dat}:
#' \itemize{
#'   \item Shared contact history
#'     (\code{common-contact-history-within-community = 1}):
#'     \code{community_id start_day stop_day contact_mode offset ignore}
#'   \item Individualised contact history (= 0):
#'     \code{person_id start_day stop_day contact_mode offset ignore}
#' }
#'
#' \strong{p2p_contact.dat}:
#' \itemize{
#'   \item Shared: \code{community_id start_day stop_day contact_mode offset ignore}
#'   \item Individualised: \code{start_day stop_day person_i person_j contact_mode offset ignore}
#' }
#'
#' \strong{impute.dat} — one row per individual requiring outcome imputation:
#' \code{person_id possible_pre_immune possible_escape
#' possible_sym possible_sym_start possible_sym_stop
#' possible_asym possible_asym_start possible_asym_stop}
#'
#' @param cfg Named list returned by \code{\link{read_config}}.  The
#'   \code{path_in} element determines the directory from which data files are
#'   read.
#' @param names_tid Character vector of covariate names for the
#'   time-independent covariate file.  Must have length equal to
#'   \code{cfg$n_time_ind_covariate}.  When \code{NULL} (the default) the
#'   covariates are named \code{x1}, \code{x2}, \ldots, \code{x\{n_tic\}}
#'   where the indices run sequentially from 1.
#' @param names_tdp Character vector of covariate names for the
#'   time-dependent covariate file.  Must have length equal to
#'   \code{cfg$n_time_dep_covariate}.  When \code{NULL} (the default) the
#'   covariates are named \code{x\{n_tic+1\}}, \ldots,
#'   \code{x\{n_tic+n_tdc\}}, continuing the sequential scheme used for the
#'   time-independent covariates.
#'
#' @return A named list with seven elements, each a data frame (or \code{NULL}
#'   if the corresponding file was absent):
#' \describe{
#'   \item{\code{pop}}{Individual-level data frame (13 columns).}
#'   \item{\code{community}}{Community-level data frame (5 columns).}
#'   \item{\code{time_ind_covariate}}{Time-independent covariate data frame
#'     (\code{NULL} when \code{n_time_ind_covariate} is 0 in the config).}
#'   \item{\code{time_dep_covariate}}{Time-dependent covariate data frame
#'     (\code{NULL} when \code{n_time_dep_covariate} is 0 in the config).}
#'   \item{\code{c2p_contact}}{Community-to-person contact data frame, or
#'     \code{NULL} if auto-generated or file absent.}
#'   \item{\code{p2p_contact}}{Person-to-person contact data frame, or
#'     \code{NULL} if auto-generated or file absent.}
#'   \item{\code{impute}}{Outcome-imputation data frame (9 columns), or
#'     \code{NULL} if the EM option is off or file absent.}
#' }
#'
#' @seealso \code{\link{gen_population}}, \code{\link{write_population}},
#'   \code{\link{read_config}}, \code{\link{simulate_single}}
#'
#' @examples
#' \dontrun{
#' cfg_file <- system.file("extdata", "CaseStudy2", "config.file",
#'                         package = "ChainBinomial")
#' cfg <- read_config(cfg_file)
#' dat <- read_population(cfg)
#'
#' # Quick summary
#' nrow(dat$pop)
#' table(dat$pop$infection)
#'
#' # First few rows of the community table
#' head(dat$community)
#' }
#'
#' @export
read_population <- function(cfg, names_tid = NULL, names_tdp = NULL) {

    path_in <- cfg$path_in

    # ------------------------------------------------------------------ #
    # Helper: read a whitespace-delimited file if it exists.              #
    # Returns NULL (invisibly) when the file is absent.                   #
    # ------------------------------------------------------------------ #
    read_if_exists <- function(fname, col_names) {
        fpath <- file.path(path_in, fname)
        if (!file.exists(fpath)) return(invisible(NULL))
        read.table(fpath, header = FALSE, col.names = col_names)
    }

    # ---- 1.  pop.dat -------------------------------------------------- #
    pop <- read_if_exists(
        "pop.dat",
        c("id", "community", "pre_immune", "infection", "symptom",
          "day_ill", "exit", "day_exit", "idx",
          "u_mode", "q_mode", "weight", "ignore")
    )

    # ---- 2.  community.dat -------------------------------------------- #
    community <- read_if_exists(
        "community.dat",
        c("community", "day_epi_start", "day_epi_stop",
          "day_last_followup", "c2p_group")
    )

    # ---- 3.  time_ind_covariate.dat ------------------------------------ #
    n_tic <- cfg$n_time_ind_covariate
    tic_names_default <- if (n_tic > 0L) paste0("x", seq_len(n_tic)) else character(0L)
    if (!is.null(names_tid)) {
        if (length(names_tid) != n_tic)
            stop("length(names_tid) = ", length(names_tid),
                 " but there are ", n_tic,
                 " time-independent covariate(s) in the config.")
        tic_col_names <- as.character(names_tid)
    } else {
        tic_col_names <- tic_names_default
    }
    time_ind_covariate <- if (n_tic > 0L) {
        df <- read_if_exists("time_ind_covariate.dat",
                             c("id", paste0("v_", seq_len(n_tic))))
        if (!is.null(df)) names(df) <- c("id", tic_col_names)
        df
    } else {
        NULL
    }

    # ---- 4.  time_dep_covariate.dat ------------------------------------ #
    n_tdc <- cfg$n_time_dep_covariate
    tdc_names_default <- if (n_tdc > 0L) paste0("x", n_tic + seq_len(n_tdc)) else character(0L)
    if (!is.null(names_tdp)) {
        if (length(names_tdp) != n_tdc)
            stop("length(names_tdp) = ", length(names_tdp),
                 " but there are ", n_tdc,
                 " time-dependent covariate(s) in the config.")
        tdc_col_names <- as.character(names_tdp)
    } else {
        tdc_col_names <- tdc_names_default
    }
    time_dep_covariate <- if (n_tdc > 0L) {
        df <- read_if_exists("time_dep_covariate.dat",
                             c("id", "day_start", "day_stop",
                               paste0("v_", seq_len(n_tdc))))
        if (!is.null(df)) names(df) <- c("id", "day_start", "day_stop", tdc_col_names)
        df
    } else {
        NULL
    }

    # ---- 5.  c2p_contact.dat ------------------------------------------ #
    # Shared history (common_contact_history_within_community == 1):
    #   community_id  start_day  stop_day  contact_mode  offset  ignore
    # Individualised history (== 0):
    #   person_id  start_day  stop_day  contact_mode  offset  ignore
    c2p_cols <- if (cfg$common_contact_history_within_community == 1L) {
        c("community_id", "start_day", "stop_day", "contact_mode", "offset", "ignore")
    } else {
        c("person_id", "start_day", "stop_day", "contact_mode", "offset", "ignore")
    }
    c2p_contact <- read_if_exists("c2p_contact.dat", c2p_cols)

    # ---- 6.  p2p_contact.dat ------------------------------------------ #
    # Shared history (== 1):
    #   community_id  start_day  stop_day  contact_mode  offset  ignore
    # Individualised history (== 0):
    #   start_day  stop_day  person_i  person_j  contact_mode  offset  ignore
    p2p_cols <- if (cfg$common_contact_history_within_community == 1L) {
        c("community_id", "start_day", "stop_day", "contact_mode", "offset", "ignore")
    } else {
        c("start_day", "stop_day", "person_i", "person_j", "contact_mode", "offset", "ignore")
    }
    p2p_contact <- read_if_exists("p2p_contact.dat", p2p_cols)

    # ---- 7.  impute.dat ----------------------------------------------- #
    impute <- read_if_exists(
        "impute.dat",
        c("person_id",
          "possible_pre_immune", "possible_escape",
          "possible_sym",  "possible_sym_start",  "possible_sym_stop",
          "possible_asym", "possible_asym_start", "possible_asym_stop")
    )

    # ------------------------------------------------------------------ #
    list(
        pop                = pop,
        community          = community,
        time_ind_covariate = time_ind_covariate,
        time_dep_covariate = time_dep_covariate,
        c2p_contact        = c2p_contact,
        p2p_contact        = p2p_contact,
        impute             = impute
    )
}
