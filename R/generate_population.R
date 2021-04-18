#' @name generate_population
#' @title Generate population of animals with genetic structure and allocated
#' in groups
#'
#' @import MASS
#'
#' @param num_replications number of replications
#' @param sires number of sires
#' @param dpsire number of dam per sires
#' @param group_size group size
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
#' grp_size = 10
#' example <- generate_population(num_replications = 10, group_size = grp_size, sires = 25, dpsire = 2,
#'             rhoG = -0.4, rhoE = -0.4, SigG.g = 4, SigG.f = 4, SigE.g = 1,
#'              SigE.f = 1)
#' example
#'
#'
#' @export
library(MASS)

# assumes same relationship matrix for all replicates!
generate_population <- function(num_replications = 1, sires = 100, dpsire = 20, rhoG = 0,
  rhoE = 0, SigG.g = 4, SigG.f = 4, SigE.g = 1, SigE.f = 1, group_size = 10, seed = 242) {

  N <- sires * dpsire
  ngroups = N/group_size
  n.family = sires   # número de famílias
  size.family =  dpsire   # tamanho da família
  ng.family=dpsire/(0.5*group_size)    # número de grupos por família
  ng.block = (ng.family*(1 + ng.family))/2   # número de grupos por bloco
  nf.block = (ng.family + 1)    # número de famílias por bloco
  n.blocks = (n.family*size.family)/(ng.block*group_size)   # número de blocos
  
  set.seed(seed)

  BV_sire_all_replicates <- data.frame(replicate = numeric(), sire_ID = numeric(), ng.off = numeric(),
                         Ag = numeric(), Af = numeric(), stringsAsFactors = FALSE)

  BV_dam_all_replicates <- data.frame(replicate = numeric(), dam_ID = numeric(), sire_ID = numeric(),
                         Ag = numeric(), Af = numeric(), stringsAsFactors = FALSE)

  offspring_all_replicates <- data.frame(replicate = numeric(), ID = numeric(), sire_ID = numeric(),
                              group = numeric(), index = numeric(), Ag = numeric(),
                              Af = numeric(), Eg = numeric(), Ef = numeric(),
                              stringsAsFactors = FALSE)

  for (replic in 1:num_replications) {

    #----------------------------------------#
    # definindo os valores genéticos #
    #----------------------------------------#
    # parents
    covG = rhoG * sqrt(SigG.g) * sqrt(SigG.f)
    sigmaG = matrix(c(SigG.g, covG, covG, SigG.f), nc = 2, byrow = T)
    auxBV_sire <- mvrnorm(n = sires, mu = c(0, 0), Sigma = sigmaG)
    auxBV_dam <- mvrnorm(n = N, mu = c(0, 0), Sigma = sigmaG)

    BV_sire <- data.frame(replicate = replic, sire_ID = 1:sires, ng.off = NA, Ag = auxBV_sire[, 1],
      Af = auxBV_sire[, 2])
    BV_dam <- data.frame(replicate = replic, dam_ID = 1:N, sire_ID = sort(rep(1:sires, dpsire)),
      Ag = auxBV_dam[, 1], Af = auxBV_dam[, 2])

  #   # offspring
    offspring <- data.frame(replicate = replic, ID = (1:N), sire_ID = sort(rep(1:sires,
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
    environ_effects <- mvrnorm(n = N, mu = c(0, 0), Sigma = sigmaE)
    offspring$Eg <- environ_effects[, 1]
    offspring$Ef <- environ_effects[, 2]
    offspring$g <- offspring$Ag + offspring$Eg
    offspring$f <- offspring$Af + offspring$Ef

    #------------------------------------------------------#
    # Assigning index cases and groups #
    #------------------------------------------------------#
    #------------------------------#
    # Random family assignment #
    #------------------------------#
    #------------------------------------#
    # Formando os blocos                #
    # Gerando os Grupos Aleatoriamente #
    #----------------------------------#
    
    for (block in 1:n.blocks) {
      
      # Identificando a amplitude dos pais
      p1 <- 5*(block-1) +1 
      p2 <- 5*block
      
      # Criando as combinações de famílias
      comb_pais <- combn(seq(p1, p2, 1), 2)
      
      # Seleção aleatoria para cada grupo
      for (jj in 1:group_size) {
        combinacao <- comb_pais[, jj]
        
        indiv1 <- sample(rownames(offspring[offspring$sire_ID == combinacao[1] & is.na(offspring$groups), ]), 5, replace = F)
        indiv2 <- sample(rownames(offspring[offspring$sire_ID == combinacao[2] & is.na(offspring$groups), ]), 5, replace = F)
        
        
        offspring[c(indiv1, indiv2), "groups"] <- paste(combinacao[1], combinacao[2], sep = "|")
        
      }
    }
    
    
    
    # index cases
    #offspring[, "index"] <- sample(c(rep(0,ngroups),rep(1,N-ngroups)),replace=F)
    # groups:
    #offspring$group <- NA
    #offspring[offspring$index == 0, ]$group <- 1:ngroups
    #offspring[offspring$index == 1, ]$group <- sample(rep(1:ngroups, group_size -       1), replace = F)
    #offspring[offspring$index == 0, 'tau'] <- 0

  #   # TODO include more than one replication in the same data frame
    new.order <- c("replicate","ID", "sire_ID", "group", "index", "Ag", "Af", "Eg",
      "Ef")
    offspring <- offspring[, new.order]

    for (i in 1:sires) BV_sire[i, "ng.off"] <- with(offspring[offspring$sire_ID ==
      i, ], dim(table(group)))


    BV_sire_all_replicates  <- rbind(BV_sire_all_replicates, BV_sire)
    BV_dam_all_replicates  <- rbind(BV_dam_all_replicates, BV_dam)
    offspring_all_replicates <- rbind(offspring_all_replicates, offspring)
    }

    # relationship matrix (A)
    relationship_matrix_aux <-list()
    relationship_matrix_aux[[1]] <- diag(sires)
    for (i in 2:(sires+1)){
      relationship_matrix_aux[[i]] <- matrix(0.25, dpsire, dpsire)
      diag(relationship_matrix_aux[[i]]) <- 1
    }
    relationship_matrix <- dlm::bdiag(relationship_matrix_aux)



  return(list(sire = BV_sire_all_replicates,
              dam = BV_dam_all_replicates,
              offspring = offspring_all_replicates,
              relationship_matrix = as.matrix(relationship_matrix)))
}
