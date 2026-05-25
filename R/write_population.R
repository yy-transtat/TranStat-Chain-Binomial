#' Write Pseudo-population to Text Files
#'
#' Writes the data frames in \code{pop_list} to text files in the format
#' expected by the TranStat C programs (\code{pop.dat},
#' \code{community.dat}, \code{time_ind_covariate.dat},
#' \code{time_dep_covariate.dat}).  Optional files are written when the
#' corresponding list element is non-\code{NULL}:
#' \code{c2p_contact.dat}, \code{p2p_contact.dat} (contact histories), and
#' \code{impute.dat} (outcome imputation, as produced by
#' \code{\link{gen_impute}}).
#'
#' @param pop_list Named list returned by \code{\link{gen_population}}.
#' @param dir Character. Directory in which to write the files. Created if it
#'   does not exist. Default is the current working directory.
#'
#' @return Called for its side effect of writing files. Returns \code{dir}
#'   invisibly.
#'
#' @examples
#' \dontrun{
#' pop <- gen_population(
#'   n_community    = 100,
#'   community_size = 5,
#'   day_epi_stop   = 14,
#'   case_ascertained = 1L
#' )
#' write_population(pop, dir = "data/example")
#' }
#'
#' @export
write_population <- function(pop_list, dir = ".") {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)

    ## ---- pop.dat ----
    ## format: %8d  %8d  %1d  %1d  %1d  %8d  %1d  %8d  %1d  %2d  %2d  %f %1d
    ## The C engine uses 99999 as the sentinel for missing day_ill / day_exit;
    ## R represents the same condition as NA, so replace before formatting.
    p <- pop_list$pop
    p$day_ill[is.na(p$day_ill)]   <- 99999L
    p$day_exit[is.na(p$day_exit)] <- 99999L
    pop_lines <- sprintf(
        "%8d  %8d  %1d  %1d  %1d  %8d  %1d  %8d  %1d  %2d  %2d  %f %1d",
        p$id, p$community, p$pre_immune, p$infection, p$symptom,
        p$day_ill, p$exit, p$day_exit, p$idx,
        p$u_mode, p$q_mode, p$weight, p$ignore
    )
    writeLines(pop_lines, file.path(dir, "pop.dat"))

    ## ---- community.dat ----
    ## format: %8d  %8d  %8d  %8d  %8d
    cm <- pop_list$community
    com_lines <- sprintf(
        "%8d  %8d  %8d  %8d  %8d",
        cm$community, cm$day_epi_start, cm$day_epi_stop,
        cm$day_last_followup, cm$c2p_group
    )
    writeLines(com_lines, file.path(dir, "community.dat"))

    ## ---- time_ind_covariate.dat (optional) ----
    ## format: %8d  %f
    if (!is.null(pop_list$time_ind_covariate)) {
        tic <- pop_list$time_ind_covariate
        tic_lines <- sprintf("%8d  %f", tic$id, tic$value)
        writeLines(tic_lines, file.path(dir, "time_ind_covariate.dat"))
    }

    ## ---- time_dep_covariate.dat (optional) ----
    ## format: %8d  %8d  %8d  %f  (id, day_start, day_stop, value)
    if (!is.null(pop_list$time_dep_covariate)) {
        tdc <- pop_list$time_dep_covariate
        tdc_lines <- sprintf(
            "%8d  %8d  %8d  %f",
            tdc$id, tdc$day_start, tdc$day_stop, tdc$value
        )
        writeLines(tdc_lines, file.path(dir, "time_dep_covariate.dat"))
    }

    ## ---- c2p_contact.dat / p2p_contact.dat (optional) ----
    ## Written only when the corresponding data frame is non-NULL.
    ## Column layout (both shared and individualised variants):
    ##   6-column c2p / shared p2p:
    ##     id  start_day  stop_day  contact_mode  offset  ignore
    ##     format: %8d  %8d  %8d  %2d  %f  %1d
    ##   7-column individualised p2p:
    ##     start_day  stop_day  person_i  person_j  contact_mode  offset  ignore
    ##     format: %8d  %8d  %8d  %8d  %2d  %f  %1d
    ## The offset column (second-to-last) is the only double; a local helper
    ## formats each column individually to handle both layouts uniformly.
    .write_contact_dat <- function(df, fname) {
        n   <- ncol(df)
        off <- n - 1L           # offset column (1-based): second-to-last
        cols <- lapply(seq_len(n), function(j) {
            if      (j == off) sprintf("%f",  df[[j]])   # offset: double
            else if (j == n)   sprintf("%1d", df[[j]])   # ignore flag
            else if (j == n - 2L && n == 7L)             # contact_mode (7-col layout)
                               sprintf("%2d", df[[j]])
            else if (j == 4L && n == 6L)                 # contact_mode (6-col layout)
                               sprintf("%2d", df[[j]])
            else               sprintf("%8d", df[[j]])   # all other integers
        })
        writeLines(do.call(paste, c(cols, list(sep = "  "))), fname)
    }
    if (!is.null(pop_list$c2p_contact))
        .write_contact_dat(pop_list$c2p_contact, file.path(dir, "c2p_contact.dat"))
    if (!is.null(pop_list$p2p_contact))
        .write_contact_dat(pop_list$p2p_contact, file.path(dir, "p2p_contact.dat"))

    ## ---- impute.dat (optional) ----
    ## Written only when pop_list$impute is non-NULL (asymptomatic infections
    ## present).  Format mirrors core.h fscanf calls:
    ##   person_id  possible_pre_immune  possible_escape
    ##   possible_sym  possible_sym_start  possible_sym_stop
    ##   possible_asym  possible_asym_start  possible_asym_stop
    if (!is.null(pop_list$impute)) {
        imp <- pop_list$impute
        imp_lines <- sprintf(
            "%8d  %1d  %1d  %1d  %8d  %8d  %1d  %8d  %8d",
            imp$person_id,
            imp$possible_pre_immune, imp$possible_escape,
            imp$possible_sym, imp$possible_sym_start, imp$possible_sym_stop,
            imp$possible_asym, imp$possible_asym_start, imp$possible_asym_stop
        )
        writeLines(imp_lines, file.path(dir, "impute.dat"))
    }

    invisible(dir)
}
