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
  vector[2] mean_A; 
  vector[2] mean_E;
  }


parameters {
  cholesky_factor_corr[2] ro_A;
  vector<lower=0>[2] sigma_s;
  cholesky_factor_corr[2] ro_E;
  vector<lower=0>[2] sigma_N;
  matrix[S,2] AS; 
  matrix[N,2] EN;
}

transformed parameters{
  matrix[2,2] Omega_A;
  matrix[2,2] Sigma_A;
  matrix[2,2] Omega_E;
  matrix[2,2] Sigma_E;
  Omega_A = multiply_lower_tri_self_transpose(ro_A);
  Sigma_A = quad_form_diag(Omega_A, sigma_s);
  Omega_E = multiply_lower_tri_self_transpose(ro_E);
  Sigma_E = quad_form_diag(Omega_E, sigma_N);
}

model {
  real lambda;
  real infect1;
  real infect2;

  target += lkj_corr_cholesky_lpdf(ro_A | 1);
  for (i in 1:2){
    target += cauchy_lpdf(sigma_s[i] | 0, 5);
  }
  target += lkj_corr_cholesky_lpdf(ro_E | 1);
  for (i in 1:2){
    target += cauchy_lpdf(sigma_N[i] | 0, 5);
  }
  
  for (i in 1:S){
    target += multi_normal_lpdf(AS[i,] | mean_A, Sigma_A); 
    }

  for (i in 1:N){
    if (n_index[i]==0 && is_last_infection[i]==0){ // contribui apenas com a infectividade
      target += normal_lpdf(EN[i,2] | mean_E[2], sqrt(Sigma_E[2,2])); 
    }
    if (n_index[i]==1 && is_last_infection[i]==0){ // contribui com infectividade e susceptibilidade
      target += multi_normal_lpdf(EN[i,] | mean_E, Sigma_E); 
    }
    if (n_index[i]==1 && is_last_infection[i]==1){ // contribui apenas com susceptibilidade
      target += normal_lpdf(EN[i,1] | mean_E[1], sqrt(Sigma_E[1,1])); 
    }
  }

  for(j in 1:N){
    if(n_index[j]==1){
      lambda = 0;
      infect1 = 0;
      infect2 = 0;
      target += AS[sire_ID[j],1] + EN[j,1];
      for (k in ((group[j]-1)*group_size+1):(j-1)){
      	infect1 = infect1+exp(AS[sire_ID[k],2] + EN[k,2]);
      	infect2 = infect2+((infection_time[j]-infection_time[k])*exp(AS[sire_ID[k],2] + EN[k,2]));
      } 
      target += log(infect1)-(exp(AS[sire_ID[j],1] + EN[j,1])*infect2);
    }
  }
}
    
