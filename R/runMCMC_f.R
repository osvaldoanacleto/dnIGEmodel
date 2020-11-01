runMCMC_f <- function(xvector) {
  if (!is.double(xvector)) {storage.mode(xvector) <- 'double'}
    .Call(c_runMCMC_f, xvector)
}
