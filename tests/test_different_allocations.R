library(devtools)
#install_github("osvaldoanacleto/dnIGEmodel", ref = "experimental_design")
library(dnIGEmodel)

grp_size <- 10

result_random <- generate_population(
    allocation_type = "random",
    num_replications =1,
    group_size = grp_size,
    sires = 100,
    dpsire = 20
  )

result_2FAM <-
  generate_population(
    allocation_type = "2FAM",
    num_replications = 1,
    group_size = grp_size,
    sires = 100,
    dpsire = 20,
  )


result_random
result_2FAM