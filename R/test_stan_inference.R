# TODO: add dependencies to package
library(rstan)
library(ggplot2)
library(compiler)


test_stan_inference <- function(num_replications,
                                group_size,
                                sires,
                                dpsire) {
  offspring <-
    generate_population(
      num_replications = num_replications,
      group_size =  group_size,
      sires = sires,
      dpsire = dpsire
    )$offspring
  
  offspring_epidemics <- transform(
    generate_epidemics(
      offspring_data_replicates =  offspring,
      group_size = group_size
    ),
    time_between_infections = ave(
      infection_time,
      replicate,
      group,
      FUN = function(x)
        c(0,  diff(x))
    )

  )
  
  estimates <-
    estimate_model_parameters_stan(
      offspring_replicates = offspring_epidemics,
      num_replications = num_replications,
      N = sires * dpsire,
      S = sires,
      ngroups = N/group_size,
      n_chains = 1,
      n_iterations = 1000,
      burnin_size = 100,
      thin_size = 10
    )
  
  # inference_results <- generate_inference_results()
  
  return(estimates)
  
}

results <- test_stan_inference(
  num_replications = 1,
  group_size = 10,
  sires = 50,
  dpsire = 10
)


# variance_chains <- as.data.frame(results)
summary(results)
