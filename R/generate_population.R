#' @name generate_population
#' @title Generate population of animals with genetic structure and allocated
#' in groups
#'
#' @import MASS
#'
#' @param sires number of sires
#' @param dpsire number of dam per sires
#' @param gr.size group size
#' @param index 1 if animal is index case, 0 otherwise
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
    auxBV_sire <- mvrnorm(n = sires, mu = c(0, 0), Sigma = sigmaG)
    auxBV_dam <- mvrnorm(n = N, mu = c(0, 0), Sigma = sigmaG)

    BV_sire <- data.frame(sire_ID = 1:sires, ng.off = NA, Ag = auxBV_sire[, 1],
      Af = auxBV_sire[, 2])
    BV_dam <- data.frame(dam_ID = 1:N, sire_ID = sort(rep(1:sires, dpsire)),
      Ag = auxBV_dam[, 1], Af = auxBV_dam[, 2])

    # offspring
    offspring <- data.frame(offspring_ID = (1:N), sire_ID = sort(rep(1:sires,
      dpsire)), Ag = rep(NA, N), Af = rep(NA, N))
    for (i in 1:sires) {
      mendelian_term <- mvrnorm(n = dpsire, mu = c(0, 0), Sigma = 0.5 * sigmaG)
      offspring[offspring$sire_ID == i, ]$Ag <- 0.5 * BV_sire$Ag[i] + 0.5 *
        BV_dam[BV_dam$sire_ID == i, ]$Ag + mendelian_term[, 1]
      offspring[offspring$sire_ID == i, ]$Af <- 0.5 * BV_sire$Af[i] + 0.5 *
        BV_dam[BV_dam$sire_ID == i, ]$Af + mendelian_term[, 2]
    }

    #------------------------------------------------------#
    # Simulating Phenotypes #
    #------------------------------------------------------#
    covE = rhoE * sqrt(SigE.g) * sqrt(SigE.f)
    sigmaE = matrix(c(SigE.g, covE, covE, SigE.f), nc = 2, byrow = T)
    Eaux <- mvrnorm(n = N, mu = c(0, 0), Sigma = sigmaE)
    offspring$Eg <- Eaux[, 1]
    offspring$Ef <- Eaux[, 2]
    offspring$g <- offspring$Ag + offspring$Eg
    offspring$f <- offspring$Af + offspring$Ef

    #------------------------------------------------------#
    # Assigning index cases and groups #
    #------------------------------------------------------#
    #------------------------------#
    # Random family assignment #
    #------------------------------#
    # index cases
    offspring[, "index"] <- sample(c(rep(1, ngroups), rep(0, N - ngroups)), replace = F)
    # groups:
    offspring$group <- NA
    offspring[offspring$index == 1, ]$group <- 1:ngroups
    offspring[offspring$index == 0, ]$group <- sample(rep(1:ngroups, gr.size -
      1), replace = F)
    # offspring[offspring$s == 0, 'tau'] <- 0

    # TODO include more than one replication in the same data frame
    new.order <- c("offspring_ID", "sire_ID", "group", "index", "Ag", "Af", "Eg",
      "Ef")
    offspring <- offspring[, new.order]

    for (i in 1:sires) BV_sire[i, "ng.off"] <- with(offspring[offspring$sire_ID ==
      i, ], dim(table(group)))
  }

  return(list(sire = BV_sire, dam = BV_dam, offspring = offspring))
}
