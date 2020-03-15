#' @name generate_population
#' @title Generate population of animals with genetic structure and allocated
#' in groups
#'
#' @param sires number of sires
#' @param dpsire number of dam per sires
#' @param gr.size group size
#' @param SigG.g = susceptibility genetic variance
#' @param SigG.f = infectivity genetic variance
#' @param SigE.g = susceptibility environmental variance
#' @param SigE.f = infectivity environmental variance
#' @param rhoG genetic correlation between susceptibility and infectivity
#' @param rhoE environmental correlation between susceptibility and infectivity
#'
#' @return data frame with
#'
#' @examples
#' function(num_replications = 1, sires = 2, dpsire = 5,
#'          rhoG = -0.4, rhoE = -0.4, SigG.g = 4, SigG.f = 4, SigE.g = 1,
#'          SigE.f = 1,
#'          gr.size = 5, seed = 242)
#'
#' @export
library(MASS)


generate_population <- function(num_replications = 1, sires = 2, dpsire = 5, rhoG = -0.4,
  rhoE = -0.4, SigG.g = 4, SigG.f = 4, SigE.g = 1, SigE.f = 1, gr.size = 5, seed = 242) {

  N <- sires * dpsire
  ngroups = N/gr.size
  set.seed(seed)

  for (jj in 1:num_replications) {

    #----------------------------------------#
    # generating BVs #
    #----------------------------------------#
    # parents
    covG = rhoG * sqrt(SigG.g) * sqrt(SigG.f)
    sigmaG = matrix(c(SigG.g, covG, covG, SigG.f), nc = 2, byrow = T)
    auxBV.s <- mvrnorm(n = sires, mu = c(0, 0), Sigma = sigmaG)
    auxBV.d <- mvrnorm(n = N, mu = c(0, 0), Sigma = sigmaG)

    BV.s <- data.frame(animal_ID = 1:sires, ng.off = NA, Ag = auxBV.s[, 1], Af = auxBV.s[,
      2])
    BV.d <- data.frame(animal_ID = 1:N, sire = sort(rep(1:sires, dpsire)), Ag = auxBV.d[,
      1], Af = auxBV.d[, 2])

    # offspring
    pop <- data.frame(animal_ID = (1:N), sire = sort(rep(1:sires, dpsire)), Ag = rep(NA,
      N), Af = rep(NA, N))
    for (i in 1:sires) {
      auxBVmend <- mvrnorm(n = dpsire, mu = c(0, 0), Sigma = 0.5 * sigmaG)
      pop[pop$sire == i, ]$Ag <- 0.5 * BV.s$Ag[i] + 0.5 * BV.d[BV.d$sire ==
        i, ]$Ag + auxBVmend[, 1]
      pop[pop$sire == i, ]$Af <- 0.5 * BV.s$Af[i] + 0.5 * BV.d[BV.d$sire ==
        i, ]$Af + auxBVmend[, 2]
    }

    #------------------------------------------------------#
    # Simulating Phenotypes #
    #------------------------------------------------------#
    covE = rhoE * sqrt(SigE.g) * sqrt(SigE.f)
    sigmaE = matrix(c(SigE.g, covE, covE, SigE.f), nc = 2, byrow = T)
    Eaux <- mvrnorm(n = N, mu = c(0, 0), Sigma = sigmaE)
    pop$Eg <- Eaux[, 1]
    pop$Ef <- Eaux[, 2]
    pop$g <- pop$Ag + pop$Eg
    pop$f <- pop$Af + pop$Ef

    #------------------------------------------------------#
    # Assigning index cases and groups #
    #------------------------------------------------------#
    #------------------------------#
    # Random family assignment #
    #------------------------------#
    # index cases
    pop[, "s"] <- sample(c(rep(0, ngroups), rep(1, N - ngroups)), replace = F)
    # groups:
    pop$group <- NA
    pop[pop$s == 0, ]$group <- 1:ngroups
    pop[pop$s == 1, ]$group <- sample(rep(1:ngroups, gr.size - 1), replace = F)
    # pop[pop$s == 0, 'tau'] <- 0

    # TODO include more than one replication in the same data frame
    new.order <- c("animal_ID", "sire", "group", "s", "Ag", "Af", "Eg", "Ef")
    pop <- pop[, new.order]

    for (i in 1:sires) BV.s[i, "ng.off"] <- with(pop[pop$sire == i, ], dim(table(group)))
  }

  return(pop)
}
