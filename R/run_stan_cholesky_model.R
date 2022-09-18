gc()
library(MASS)
source("R/generate_epidemics.R")
# Num_replications = 20
RhoG <- -0.4
RhoE <- -0.4
SigGG <- 4
SigGF <- 4
SigEG <- 1
SigEF <- 1
HG <- SigGG / (SigGG + SigEG)
HF <- HG
CovG <- RhoG * sqrt(SigGG) * sqrt(SigGF)
CovG

# Sires = 99
# Dpsire = 20
# Group_size = 20
Num_replications <- 1
Sires <- 2
Dpsire <- 10
Group_size <- 10
N <- Sires * Dpsire
ngroups <- N / Group_size
# ==================================================================================#
#--------------------------EPIDEMIA-2FAM
# ==================================================================================#
source("R/generate_population.R")

pop_2FAM <-
  generate_population(
    num_replications = Num_replications, sires = Sires, dpsire = Dpsire,
    rhoG = RhoG, rhoE = RhoE,
    SigG.g = SigGG, SigG.f = SigGF, SigE.g = SigEG, SigE.f = SigEF,
    group_size = Group_size, seed = 0702, allocation_type = "random"
  )
head(pop_2FAM)


epi_2FAM <- generate_epidemics(pop_2FAM$offspring, group_size = Group_size, seed = 0242)
head(epi_2FAM)


library(rstan)
library(compiler)
options(buildtools.check = function(action) TRUE)
set.seed(2015)


modelocholes <- "R/sire_model_cholesky.stan"
enableJIT(3)
# for (ii in 1:Num_replications) {
for (ii in 1:1) {
  # ii<-1
  fit <- stan(
    model_code = modelocholes,
    data = list(
      N = N,
      S = Sires,
      sire_ID = subset(epi_2FAM, epi_2FAM$replicate == ii)$sire_ID,
      ngroups = ngroups,
      group = subset(epi_2FAM, epi_2FAM$replicate == ii)$group,
      group_size = Group_size,
      n_index = subset(epi_2FAM, epi_2FAM$replicate == ii)$index,
      infection_time = subset(epi_2FAM, epi_2FAM$replicate == ii)$infection_time,
      is_last_infection = subset(epi_2FAM, epi_2FAM$replicate == ii)$is_last_infection,
      mean_A = c(0, 0),
      mean_E = c(0, 0)
    ),
    chains = 2,
    iter = 11000,
    warmup = 1000,
    thin = 5
  )
  samp <- rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE, include = TRUE)
  cadeia1 <- round(samp[, 1, ], 3)
  cadeia2 <- round(samp[, 2, ], 3)
  write.table(cadeia1, file = paste("~/Simulacao/cadeia1_rep", eval(ii), ".csv", sep = ""),
              row.names = FALSE, append = FALSE)
  write.table(cadeia2, file = paste("~/Simulacao/cadeia2_rep", eval(ii), ".csv", sep = ""),
              row.names = FALSE, append = FALSE)
}
