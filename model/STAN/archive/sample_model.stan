data{
  int<lower=1> N;                // number of rows
  int<lower=1> P;                // number of cols
  matrix[N,P] Y;                 // data matrix of order [N,P]
  int<lower=2> K;             // number of latent dimensions 
  int<lower=1> K_choose_2;
}
parameters{    
  vector[K_choose_2] L_lower;   // lower diagonal elements of L
  vector<lower=0>[K] L_diag;   // lower diagonal elements of L
  vector<lower=0>[P] psi;         // vector of variances
  real<lower=0>   mu_psi;
  real<lower=0>  sigma_psi;
  real mu_lt;
  real<lower=0>  sigma_lt;
}
transformed parameters{
  cholesky_factor_cov[K] L;
  cov_matrix[P] Q;   //Covariance mat
  
  for (k in 1:K) {
    L[k, k] = 1;
  }
  {
    int i;
    for (m in 2:K) {
      for (n in 1:(m - 1)) {
        L[m, n] = L_lower[i];
        L[n, m] = 0;
        i += 1;
      }
    }
  }
  Q = L*L'+diag_matrix(psi); 
}
model {
// the hyperpriors 
   mu_psi ~ cauchy(0, 1);
   sigma_psi ~ cauchy(0,1);
   mu_lt ~ cauchy(0, 1);
   sigma_lt ~ cauchy(0,1);
// the priors 
  L_diag ~ cauchy(0,3);
  L_lower ~ cauchy(mu_lt,sigma_lt);
  psi ~ cauchy(mu_psi,sigma_psi);
//The likelihood
for( j in 1:N) {
  Y[j] ~ multi_normal(rep_vector(0.0,P),Q); 
  }
}

