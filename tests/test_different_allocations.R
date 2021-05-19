<<<<<<< HEAD
=======
library(devtools)
#install_github("osvaldoanacleto/dnIGEmodel", ref = "experimental_design")
library(dnIGEmodel)

>>>>>>> 1bd3ee954e48b4bb002a45ecbd2fbdf59fc0141c
grp_size <- 10

result_random <- generate_population(
    allocation_type = "random",
    num_replications =1,
    group_size = grp_size,
<<<<<<< HEAD
    sires = 25,
    dpsire = 2
=======
    sires = 100,
    dpsire = 20
>>>>>>> 1bd3ee954e48b4bb002a45ecbd2fbdf59fc0141c
  )

result_2FAM <-
  generate_population(
    allocation_type = "2FAM",
    num_replications = 1,
    group_size = grp_size,
<<<<<<< HEAD
    sires = 25,
    dpsire = 2,
  )
=======
    sires = 100,
    dpsire = 20,
  )


result_random
result_2FAM
>>>>>>> 1bd3ee954e48b4bb002a45ecbd2fbdf59fc0141c
