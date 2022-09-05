# from https://onunicornsandgenes.blog/2019/10/13/using-r-animal-model-with-simulated-data/

library(AlphaSimR)
 

founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#Set simulation parameters
SP  <-  SimParam$new(founderPop)
SP$setGender  <- "yes_sys"

## Founder population
FOUNDERPOP <- runMacs(nInd = 100,
                      nChr = 20,
                      inbred = FALSE,
                      species = "GENERIC")
 
## Simulation parameters 
SIMPARAM <- SimParam$new(FOUNDERPOP)
SIMPARAM$addTraitA(nQtlPerChr = 100,
                   mean = 100,
                   var = 10)
SIMPARAM$setGender("yes_sys")
SIMPARAM$setVarE(h2 = 0.3)
  
## Random mating for 9 more generations
generations <- vector(mode = "list", length = 10) 
generations[[1]] <- newPop(FOUNDERPOP,
                           simParam = SIMPARAM)
 
 
for (gen in 2:10) {
 
    generations[[gen]] <- randCross(generations[[gen - 1]],
                                    nCrosses = 10,
                                    nProgeny = 10,
                                    simParam = SIMPARAM)
 
}
 
## Put them all together
combined <- Reduce(c, generations)
 
 
## Extract phentoypes
pheno <- data.frame(animal = combined@id,
                    pheno = combined@pheno[,1])
 
## Extract pedigree
ped <- data.frame(id = combined@id,
                  dam = combined@mother,
                  sire =combined@father)
ped$dam[ped$dam == 0] <- NA
ped$sire[ped$sire == 0] <- NA
 
## Write out the files
write.csv(pheno,
          file = "R/test_stan_animal_model/sim_pheno.csv",
          row.names = FALSE,
          quote = FALSE)
 
write.csv(ped,
          file = "R/test_stan_animal_model/sim_ped.csv",
          row.names = FALSE,
          quote = FALSE)

