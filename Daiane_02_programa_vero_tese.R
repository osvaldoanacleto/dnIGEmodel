modelot<-"
data {
  int<lower=0> N;  
  int<lower=0> S; 
  int<lower=1, upper=S> sire_ID[N]; 
  int<lower=0> ngroups; 
  int<lower=1, upper=ngroups> group[N]; 
  int<lower=0> group_size; 
  vector[N] n_index;
  vector[N] infection_time;
  vector[N] is_last_infection;
  real mean_Ag; 
  real mean_Af; 
  }

parameters {
  real<lower=0> sigma_Sg; 
  real<lower=0> sigma_Eg; 
  real<lower=0> sigma_Sf; 
  real<lower=0> sigma_Ef; 
  vector[S] Ag; 
  vector[N] Eg;
  vector[S] Af; 
  vector[N] Ef;
}

transformed parameters{
  real<lower=0> precision_Sg;
  real<lower=0> precision_Eg;
  real<lower=0> precision_Sf;
  real<lower=0> precision_Ef;
  precision_Sg = 1/sigma_Sg;
  precision_Eg = 1/sigma_Eg;
  precision_Sf = 1/sigma_Sf;
  precision_Ef = 1/sigma_Ef;
}

model {
  real lambda;
  real infect1;
  real infect2;

  target += gamma_lpdf(precision_Sg | 0.01, 0.01); 
  target += gamma_lpdf(precision_Eg | 0.01, 0.01); 
  target += gamma_lpdf(precision_Sf | 0.01, 0.01); 
  target += gamma_lpdf(precision_Ef | 0.01, 0.01); 

  for (i in 1:S){
    target += normal_lpdf(Ag[i] | mean_Ag, sqrt(sigma_Sg)); 
    target += normal_lpdf(Af[i] | mean_Af, sqrt(sigma_Sf)); 
    }

  for (i in 1:N){
    if (n_index[i]==1){
      target += normal_lpdf(Eg[i] | 0, sqrt(sigma_Eg)); 
    }
    if (is_last_infection[i]==0){   // ASSUMES PREVALENCE 1 
      target += normal_lpdf(Ef[i] | 0, sqrt(sigma_Ef)); 
    }
  }

  for(j in 1:N){
    if(n_index[j]==1){
      lambda = 0;
      infect1 = 0;
      infect2 = 0;
      target += Ag[sire_ID[j]] + Eg[j];
      for (k in ((group[j]-1)*group_size+1):(j-1)){
      	infect1 = infect1+exp(Af[sire_ID[k]] + Ef[k]);
      	infect2 = infect2+((infection_time[j]-infection_time[k])*exp(Af[sire_ID[k]] + Ef[k]));
      } 
      target += log(infect1)-(exp(Ag[sire_ID[j]] + Eg[j])*infect2);
    }
  }
}      
"
Num_replications = 10
Sires = 100
Dpsire = 24
Group_size = 12
RhoG = 0
RhoE = 0
SigGG = 4
SigEG = 1
SigGF = 4
SigEF = 1
N<-Sires*Dpsire
ngroups = N/Group_size
pop_random<-
  generate_population(num_replications = Num_replications, sires = Sires, dpsire = Dpsire, rhoG = RhoG,
                      rhoE = RhoE, SigG.g = SigGG, SigG.f = SigGF, SigE.g = SigEG, SigE.f = SigEF, 
                      group_size = Group_size, seed = 0702,
                      allocation_type = "random")
epi_random<-generate_epidemics(pop_random$offspring, group_size=Group_size, seed = 0242)
head(epi_random)


library(rstan)
library(compiler)      
enableJIT(3)
ii<-1

fit<-list()

Num_replications <- 20

for (ii in ii:Num_replications){
fit[[ii]]<-
  rstan::extract(
  stan( 
    model_code = modelot,
    data = list(
      N=N,
      S=Sires,
      sire_ID=subset(epi_random,epi_random$replicate==ii)$sire_ID,
      ngroups=ngroups,
      group=subset(epi_random,epi_random$replicate==ii)$group,
      group_size=Group_size,
      n_index=subset(epi_random,epi_random$replicate==ii)$index,
      infection_time=subset(epi_random,epi_random$replicate==ii)$infection_time,
      is_last_infection=subset(epi_random,epi_random$replicate==ii)$is_last_infection,
      mean_Ag=0,
      mean_Af=0),
    chains = 1,
    iter = 1000,
    warmup = 500,
    thin = 5
      )
  ,  permuted = FALSE, inc_warmup = FALSE, include = TRUE)
}




CI=





Sg <- Sf <- Eg <- Ef  <- data.frame(rep = 1:Num_replications, median = NA, IC_inf = NA, IC_sup = NA)


for ( i in 1:Num_replications){
  
  Sg[i,2:4] <- quantile(fit[[i]][,1,1], c(0.50,.05, 0.95))
  Sf[i,2:4] <- quantile(fit[[i]][,1,3], c(0.50,.05, 0.95))
  Eg[i,2:4] <- quantile(fit[[i]][,1,2], c(0.50,.05, 0.95))
  Ef[i,2:4] <- quantile(fit[[i]][,1,4], c(0.50,.05, 0.95))
    
}


Sg
Sf
Eg
Ef


# 
# #var do ag
# plot(fit[[1]][,1,1],type='l')
# # summary(fit[[1]][,1,1])
# quantile(fit[[1]][,1,1], c(0.50,.05, 0.95))
# #var do eg
# # plot(fit[[1]][,1,2],type'l')
# # summary(fit[[1]][,1,2])
# quantile(fit[[1]][,1,2], c(0.50,.05, 0.95))
# #var do af
# # plot(fit[[1]][,1,3],type='l')
# # summary(fit[[1]][,1,3])
# quantile(fit[[1]][,1,3], c(0.50,.05, 0.95))
# #var do ef
# # plot(fit[[1]][,1,4],type='l')
# # summary(fit[[1]][,1,4])
# quantile(fit[[1]][,1,4], c(0.50,.05, 0.95))


# pais_ag<-apply(samp[,1,5:(Sires+4)],2,mean)
# cor(pais_ag,pop_random$sire$Ag)
# plot(pais_ag,pop_random$sire$Ag)
# 
# pais_af<-apply(samp[,1,(Sires+N+5):(2*Sires+N+4)],2,mean)
# cor(pais_af,pop_random$sire$Af)
# plot(pais_af,pop_random$sire$Af,xlim=c(-5,5))

