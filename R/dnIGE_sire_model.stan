data {
  int<lower=0> N;  
  int<lower=0> S; 
  int<lower=1, upper=S> sire_ID[N]; 
  int<lower=0> ngroups; 
  int<lower=1, upper=ngroups> group[N]; 
  int<lower=0> group_size; 
  vector[N] n_index;
  vector[N] infection_time;
  vector[N] ID;
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
  real infectivity_term;
    
  target += gamma_lpdf(precision_Sg | 0.001, 0.001); 
  target += gamma_lpdf(precision_Eg | 0.001, 0.001); 
  target += gamma_lpdf(precision_Sf | 0.001, 0.001); 
  target += gamma_lpdf(precision_Ef | 0.001, 0.001); 
  
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
      target += (Ag[sire_ID[j]] + Eg[j]) ;
      infectivity_term = 0;
       for(k in 1:N){
        if( (group[j] == group[k]) && 
            (ID[k] != ID[j]) &&
            (infection_time[k] < infection_time[j]))
          {
          target +=  exp(Af[sire_ID[k]] + Ef[k]);
          infectivity_term = infectivity_term + 
                    ((infection_time[j] - infection_time[k]) *
                    exp(Af[sire_ID[k]] + Ef[k]));
          }
      }
        target += - exp(Ag[sire_ID[j]] + Eg[j]) * infectivity_term;

    }
  }
}


