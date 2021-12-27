stan_model_file <-
  "/Users/osvaldo/Documents/repos/dnIGEmodel/R/dnIGE_sire_model.stan"
# option  + ß: <- !
enableJIT(3)
estimate_model_parameters_stan <- function(offspring_replicates,
                                           num_replications,
                                           N,
                                           S,
                                           ngroups,
                                           n_chains,
                                           n_iterations,
                                           burnin_size,
                                           thin_size) {
  for (ii in 1:num_replications) {
    fit <- stan(
      file = stan_model_file,
      data = list(
        N = N,
        ngroups = ngroups,
        group_size = group_size,
        S = S,
        mean_Ag = 0,
        mean_Af = 0, 
        n_index = subset(
          offspring_replicates,
          offspring_replicates$replicate == ii
        )$index,
        sire_ID = subset(
          offspring_replicates,
          offspring_replicates$replicate == ii
        )$sire_ID,
        group = subset(
          offspring_replicates,
          offspring_replicates$replicate == ii
        )$group,
        is_last_infection = subset(
          offspring_replicates,
          offspring_replicates$replicate == ii
        )$is_last_infection,
        infection_time = subset(
          offspring_replicates,
          offspring_replicates$replicate == ii
        )$infection_time,
       ID = subset(
          offspring_replicates,
          offspring_replicates$replicate == ii
        )$ID
        
      ),
      chains = n_chains,
      iter = n_iterations,
      warmup = burnin_size,
      thin = thin_size
    )
    samp <- rstan::extract(fit,
                           permuted = FALSE,
                           inc_warmup = FALSE,
                           include = TRUE)
    cadeia1 <- round(samp[, 1,], 3)
    #cadeia2<-round(samp[,2,],3)
    # write.table(cadeia1,file=paste("~/Simulacao/cadeia1_rep",eval(ii),".csv",sep=""), row.names=FALSE, append=FALSE)
    # write.table(cadeia2,file=paste("~/Simulacao/cadeia2_rep",eval(ii),".csv",sep=""), row.names=FALSE, append=FALSE)
  }
  return(samp[, 1, 1:4])
}
