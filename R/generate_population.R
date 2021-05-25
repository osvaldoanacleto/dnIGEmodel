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
#library(devtools)
#install_github("osvaldoanacleto/dnIGEmodel", ref = "experimental_design")
#library(dnIGEmodel)

# Com alterações
generate_population <- function(num_replications = 1, sires = 40, dpsire = 150, rhoG = 0,
                                rhoE = 0, SigG.g = 4, SigG.f = 4, SigE.g = 1, SigE.f = 1, group_size = 60, seed = 242,
                                allocation_type = "random"){
  N <- sires * dpsire
  ngroups = N/group_size
  n.family = sires 
  size.family =  dpsire   
  set.seed(seed)
  
  BV_sire_all_replicates <- data.frame(replicate = numeric(), sire_ID = numeric(), ng.off = numeric(),
                                       Ag = numeric(), Af = numeric(), stringsAsFactors = FALSE)
  
  BV_dam_all_replicates <- data.frame(replicate = numeric(), dam_ID = numeric(), sire_ID = numeric(),
                                      Ag = numeric(), Af = numeric(), stringsAsFactors = FALSE)
  
  offspring_all_replicates <- data.frame(replicate = numeric(), ID = numeric(), sire_ID = numeric(),
                                         group = numeric(), index = numeric(), Ag = numeric(),
                                         Af = numeric(), Eg = numeric(), Ef = numeric(),
                                         stringsAsFactors = FALSE)
  
  #offspring_all_replicates <- data.frame(replicate = numeric(), ID = numeric(), sire_ID = numeric(),
  #                                       group = character(), index = numeric(), Ag = numeric(),
  #                                       Af = numeric(), Eg = numeric(), Ef = numeric(),
  #                                       stringsAsFactors = FALSE)
  
  for (replic in 1:num_replications){
    
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
    
    # offspring
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
    
#========================================================#
#-------------------------Random
#========================================================#
if (allocation_type == "random"){
  #------------------------------------------------------#
  # Assigning index cases and groups #
  #------------------------------------------------------#
  # index cases ngroups 
    offspring[, "index"] <- sample(c(rep(0,ngroups),rep(1,N-ngroups)),replace=F)
  # groups:
      offspring$group <- NA
      offspring[offspring$index == 0, ]$group <- 1:ngroups
      offspring[offspring$index == 1, ]$group <- sample(rep(1:ngroups, group_size -
                                                              1), replace = F)
}
#========================================================#
#-------------------------2FAM 
#========================================================#
if (allocation_type == "2FAM"){
  N<-(n.family-menor.family)*size.maior.family+menor.family*size.menor.family
ng.family=dpsire/(0.5*group_size)
ng.block = (ng.family*(1 + ng.family))/2  
nf.block = (ng.family + 1)   
n.blocks = (ngroups*group_size)/(ng.block*group_size)
#n0.family= n.blocks*nf.block   ### Número total de Fam (nº quebrado)
#n.family= ceiling(n0.family)    ### Número total de Fam (inteiro)
#--------- Maiores Famílias
maior.family=floor(n.blocks)*nf.block 
size.maior.family=size.family
#---------menores famílias
menor.family=(n.family)-(maior.family) 
size.menor.family=((ngroups*group_size)-(n.family-menor.family)*size.family)/(menor.family)
#----Tamanho da amostra
#N2<-(n.family-menor.family)*size.maior.family+menor.family*size.menor.family
#----------------#
# VALORES '2FAM' #
#----------------#
#N = 5760
#n.family = 40
#size.family = 150 # size family
#group_size = 60
#ngroups = 96
#ng.family= 5
#n.blocks = 6 # 6 completos e o 7 incompleto
#nf.block = 6
#ng.block = 15
#--------- Maiores Famílias
maior.family=floor(n.blocks)*nf.block 
size.maior.family=size.family
#---------menores famílias
menor.family=(n.family)-(maior.family) 
size.menor.family=((N)-(n.family-menor.family)*size.family)/(menor.family)
5760/150
#-------------------#
# Gerando os 'Pais' #
#-------------------#
    
pai <- rep(NA, N)
#n_pai <- N/n.family
    
k<-1
  for (j in 1:(n.family-menor.family)) {
    pai[seq(k,k+(size.family-1),1)] <- rep(j, size.family)
    k <- k+size.family
    }
  for (j in (n.family-3):n.family) {
      pai[seq(k,k+(size.menor.family-1),1)] <- rep(j, size.menor.family)
      k <- k+size.menor.family
    }
    
    length(pai) 
    
#-------------------------------------------#
# Gerando as 'Maes' para cada um dos 'Pais' #
#-------------------------------------------#
    
mae <- rep(NA, N)
    
    i <-1
    for(p in 1:n.family){
      if (p <= (n.family-menor.family)){
        for (m in 1:size.family) {
          mae[i] <- paste(paste("mae", m, sep = "_"), p, sep = "_pai_")
          i <- i+1
        }
      }else{
        for (m in 1:size.menor.family) {
          mae[i] <- paste(paste("mae", m, sep = "_"), p, sep = "_pai_")
          i <- i+1
        }
      }
    }
    
#mae
    
pais <- as.data.frame(cbind(pai, mae))
    
    #head(pais)
    #dim(pais)
#----------------------------------#
# Gerando um filho para cada 'Mae' # 
# Meio-irmÃ£os                      #
#----------------------------------#
    
filhos <- rep(NA, N)
    
nf <- 1
for (f in 1:(N)) {
   filhos[f] <- paste(paste("filho", nf, sep = "_"), mae[f], sep = "_")
      
      if (pais[f, 1] <=  n.family- menor.family){
      if (nf == size.family){
          nf <- 1
        }else{
          nf <- nf+1
        }
      }else{
        if (nf == size.menor.family){
          nf <- 1
        }else{
          nf <- nf+1
        }
      }
    }
    
  #filhos
    
    
#familia <- as.data.frame(cbind(pais, filhos))
#head(familia)
#dim(familia)
    
familia$groups <- rep(NA, N)
    
#----------------------------------#
# FORMAR BLOCKS                    #
# Gerando os Grupos AleatÃ³riamente #
# apÃ³s definidos os blocos         #
#----------------------------------#
    
  n.blocks= N/(ng.block*group_size)
  n.blocks
    
  for (block in 1:ceiling(n.blocks)) {
      
  # Identificando o range dos pais
  p1 <- nf.block*(block-1) +1 
  p2 <- nf.block*block
      
  fam <- as.numeric(familia[familia$pai %in% seq(p1, p2), "pai"])
  max_p <- max(fam)
  min_p <- min(fam)
      
  # Criando combinaÃ§Ãµes de familias
  comb_pais <- combn(seq(min_p, max_p, 1), 2)
      
  # SeleÃ§Ã£o aleatoria para cada grupo
      for (jj in 1:ncol(comb_pais)) {
        combinacao <- comb_pais[, jj]
        
        indiv1 <- sample(rownames(familia[familia$pai == combinacao[1] & is.na(familia$groups), ]), 0.5*group_size, replace = F)
        indiv2 <- sample(rownames(familia[familia$pai == combinacao[2] & is.na(familia$groups), ]), 0.5*group_size, replace = F)
        
        
        familia[c(indiv1, indiv2), "groups"] <- paste(combinacao[1], combinacao[2], sep = "|")
        
      }
  }  
  
  offspring$group <- familia$groups
  
  offspring$index <- rep(NA, sires*dpsire)
  
  for (grp in unique(offspring$group)) {
    offspring[offspring$group == grp,]$index <- sample(c(0,rep(1,group_size-1)),replace = F)
  }   
  
  offspring[order(offspring$groups),] 
  }
#------------------------------------------------------------------#
#-------------------------3FAM 
#------------------------------------------------------------------#   
if (allocation_type == "2FAM"){   
  
#ng.family=size.family/((1/3)*group_size)    
#ng.block = ng.family*(ng.family-1) 
#nf.block = 2*ng.family + 1 
#n.blocks = (ngroups*group_size)/(ng.block*group_size)
#---------Total de Famílias  
#n0.family= n.blocks*nf.block   ### Número total de Fam (nº quebrado)
#n.family= ceiling(n0.family)    ### Número total de Fam (inteiro)
#--------- Maiores Famílias
#maior.family=floor(n.blocks)*nf.block 
#size.maior.family=size.family
#---------menores famílias
#menor.family=(n.family)-(maior.family) 
#size.menor.family=((ngroups*group_size)-(n.family-menor.family)*size.family)/(menor.family)
#----Tamanho da amostra
#N2<-(n.family-menor.family)*size.maior.family+menor.family*size.menor.family    
}
#==============================================================    
#offspring[offspring$index == 0, 'tau'] <- 0
# TODO include more than one replication in the same data frame
new.order <- c("replicate","ID", "sire_ID", "group", "index", "Ag", "Af", "Eg", "Ef")
offspring <- offspring[, new.order]
    
for (i in 1:sires) BV_sire[i, "ng.off"] <- with(offspring[offspring$sire_ID ==i, ], dim(table(group)))
    
    
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
              relationship_matrix = as.matrix(relationship_matrix))) #,data_set_familia = familia))
  result_random <- generate_population(allocation_type = "random")
  result_2FAM <- generate_population(allocation_type = "2FAM")
  #result_3FAM <- generate_population(allocation_type = "3FAM")
}

#levels(as.factor(result_random$offspring$replicate))
#head(result_random$offspring)

result_random
result_2FAM
#result_3FAM




