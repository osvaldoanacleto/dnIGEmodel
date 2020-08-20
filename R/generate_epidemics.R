#' @name generate_epidemics
#' @title Generate epidemics using the Gillespie algorithm
#' @title assumptions:
#' @title * animals with genetic structure and allocated in groups
#' @title * SI model
#' @title * all individuals are infected during the epidemic (prevalence = 1)
#'
#' @import MASS
#'
#' @param offspring_data dataframe created by generate_population function
#' @param group_size group size
#'
#' @return return offspring_data dataframe with column containing time of
#' @return infection
#'
#' @examples
#'grp_size = 10
#'df <- generate_population(group_size = grp_size, sires = 25, dpsire = 2)$offspring
#'df_epi <- generate_epidemics(df, group_size = grp_size)
#'
#' @export
library(MASS)

generate_epidemics <- function(offspring_data, group_size, seed = 242) {
  set.seed(seed)
  number_of_groups <- length(table(offspring_data$group))
  offspring_data$is_infected <- 1 - offspring_data$index
  #offspring_data[offspring_data$index == 1, "infection_time"] <- -1
  offspring_data[, "infection_time"] <- NA
  offspring_data[offspring_data$index == 0, "infection_time"] <- 0
  #offspring_data$is_last_infection <- 0

  lambda <- list()
  for (k in 1:number_of_groups) lambda[[k]] <- rep(NA, nr = group_size)

  # Gillespie algorithm (runs  independently by each group) --------------------
  for (i in 1:number_of_groups) {
    next_infection_time <- 0
    data_group_i <- offspring_data[offspring_data$group == i, ]
    while (sum(data_group_i$is_infected)/group_size < 1) {
      for (j in 1:group_size) {
        lambda[[i]][j] = 0
        if (data_group_i$is_infected[j] == 0) {
          lambda[[i]][j] = exp(data_group_i$Ag[j] + data_group_i$Eg[j]) *
                           sum(exp(data_group_i[data_group_i$is_infected == 1, ]$Af +
                                   data_group_i[data_group_i$is_infected == 1, ]$Ef))
        }
      }
      next_infection_time <- next_infection_time + rexp(1, rate = sum(lambda[[i]]))
      ID_next_inf <- sample(data_group_i$ID, size = 1, prob = lambda[[i]]/sum(lambda[[i]]))
      data_group_i[data_group_i$ID == ID_next_inf, "is_infected"] <- 1
      data_group_i[data_group_i$ID == ID_next_inf, "infection_time"] <- next_infection_time
    }
    #data_group_i[data_group_i$infection_time == max(data_group_i$infection_time), "is_last_infection"] <- 1
    offspring_data[offspring_data$group == i, ] <- data_group_i
  }

  return(offspring_data[order(offspring_data$group, offspring_data$infection_time ),
                        !(names(offspring_data) %in% "is_infected")])
}
