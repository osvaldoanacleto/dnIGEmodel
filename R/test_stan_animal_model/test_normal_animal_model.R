library(AGHmatrix)
library(rstan)
 
ped[is.na(ped)] <- 0
A  <- Amatrix(ped)

Z0  <- diag(1000)
L <- t(chol(A))
Z  <- Z0 %*% L
X <- model.matrix(~1, pheno)
 


# real h2;
# h2 = sigma_U / (sigma_U + sigma_E)


pheno$scaled_pheno <- as.vector(scale(pheno$pheno))
 
model_stan <- stan(file = 'R/test_stan_animal_model/animal_model.stan',
                   data = list(Y = pheno$scaled_pheno,
                               X = X,
                               A = A,
                               Z = Z0,
                               J = 1,
                               K = 1000,
                               N = 1000),
                               chains = 1)
 
est_h2_stan <- summary(model_stan, pars = "h2")$summary


