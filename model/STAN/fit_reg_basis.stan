data {
  int num_data;
  int num_basis;
  int num_years;
  vector[num_data] Y;
  vector[num_years] X;
  matrix[num_basis, num_years] B;
  int matchmethod[num_data] ; // method indexing
  int matchyears[num_data]; // year indexing
  int matchsubnat[num_data]; // subnat indexing 
  int M_count; // Number of methods
  int P_count; // Number of provinces
  int S_count; // Number of sectors
  int<lower=0> v ; 
  int<lower=0> s_sq ; 
}
parameters {
  vector[M_count] a0[P_count];
  vector<lower=0>[M_count] sigma_y;
  vector<lower=0>[M_count] tau;
  real<lower=0> sigma_tau;  // noise std
  vector<lower=0>[M_count] cstar; # slab scale
  vector<lower=0>[num_basis] lambda[P_count, M_count]; // local shrinkage parameter
}
transformed parameters {
  row_vector[num_basis] a[P_count, M_count];
  vector[num_years] Y_hat[P_count, M_count];
  matrix<lower=0, upper=1>[S_count, num_years] P[P_count, M_count]; 
  vector[num_basis] lambda_tilde[P_count, M_count];
  row_vector[num_basis] a_raw[P_count, M_count];
  
  for(p in 1:P_count) {
    for(m in 1:M_count) { 
      for(k in 1:num_basis){
        lambda_tilde[p,m,k]  = sqrt((cstar[m]*pow(lambda[p,m,k],2))/ ( cstar[m] + pow(sigma_tau,2)*pow(lambda[p,m,k], 2)));
        a_raw[p,m,k] = sigma_tau*lambda_tilde[p,m,k];
      }
      a[p,m,1] = a_raw[p,m,1]*tau[m];
      for (i in 2:(num_basis-1)) {
        a[p,m,i] = a[p,m, i-1] + a_raw[p,m,i]*tau[m]; // penalising splines
      }
      a[p,m,num_basis] = -sum(a[p,m,1:(num_basis-1)]); # hard sum to 0 constraint
      Y_hat[p,m,1:num_years] = a0[p,m]*X + to_vector(a[p,m, 1:num_basis]*B);
      P[p, m, 1, 1:num_years] = to_row_vector(inv_logit(Y_hat[p,m,1:num_years])) ;
      P[p, m, 2, 1:num_years] = to_row_vector(1 - P[p, m, 1, 1:num_years]) ;
    }
  }
}
model {
  for(p in 1:P_count) {
    for(m in 1:M_count){ 
      lambda[p,m] ~ cauchy(0,1);
    }
  }
  tau ~ cauchy(0, 1);
  sigma_y ~ cauchy(0, 1);
  cstar ~ inv_gamma(0.5*v, 0.5*v*s_sq);
  sigma_tau ~ cauchy(0,0.01);
  
  // Likelihood
  for (k in 1:num_data) {
    Y[k] ~ normal(Y_hat[matchsubnat[k], matchmethod[k], matchyears[k]], sigma_y[matchmethod[k]]);
  }
}

