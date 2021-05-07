
#' @name runMCMC
#' @title Estimate the parameters of the dnIGE model using fortran-based MCMC algorithm
#' 
#'
#' @import MASS
#'
#' @param num_replications number of replications
#' @param num_replications number of replications
#'
#' @return todo
#'
#' @examples
#' example <- runMCMC_f
#'
#'
#' @export

# runMCMC_f <- function(N) {
#   if (!is.double(N)) {storage.mode(N) <- 'double'}
#     .Call(c_runMCMC_f, N)
# }

runMCMC <- function(cte, xvector) {
  if (!is.double(xvector)) {storage.mode(xvector) <- 'double'}
  if (!is.integer(cte)) {storage.mode(cte) <- 'integer'}

    .Call(c_runMCMC_f, cte, xvector)
}
