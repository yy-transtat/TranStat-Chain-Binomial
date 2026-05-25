#' Estimate Chain Binomial Transmission Parameters
#'
#' Reads a \code{config.file} and the associated data files it points to,
#' runs the TranStat maximum likelihood (or EM) estimation algorithm, and
#' writes results to the output directory specified in \code{config.file}.
#'
#' The \code{config.file} controls all model settings: transmission parameters,
#' covariate structure, incubation and infectious period distributions,
#' optimization method, and file paths for input data and output. See the
#' TranStat documentation for the full format.
#'
#' @param config_file Character. Full path to the \code{config.file} that
#'   specifies model settings and file paths.
#' @param serial_number Integer. Identifier for this run, appended to some
#'   output file names when running multiple replicates. Default 1L.
#' @param seed Integer. Random seed for Monte Carlo EM sampling procedures.
#'   Default 12345678L.
#'
#' @return Called for its side effect of writing results to the output
#'   directory defined in \code{config.file}. Returns \code{NULL} invisibly.
#'
#' @examples
#' \dontrun{
#' transtat(
#'   config_file   = "/path/to/study/config.file",
#'   serial_number = 1L,
#'   seed          = 12345678L
#' )
#' }
#'
#' @export
transtat <- function(config_file, serial_number = 1L, seed = 12345678L) {
    invisible(.Call("r_transtat",
                    as.character(config_file),
                    as.integer(serial_number),
                    as.integer(seed)))
}
