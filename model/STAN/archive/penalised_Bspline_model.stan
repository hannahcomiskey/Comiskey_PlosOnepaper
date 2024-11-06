data{  // The input data
  int<lower=1> n_years; // Number of years
  int<lower=1> n_obs; // Number of observations
  int<lower=1> H; // Number of knots
  int<lower=1> K; // Number of spline cofficients (H+1)
  int<lower=1> P_count; // Number of provinces
  int<lower=1> C_count; // Number of countries
  int<lower=1> M_count; // Number of methods
  int<lower=1> S_count; // Number of sectors
  matrix[n_years, K] Bik; // Basis functions
  int matchcountry[P_count]; // country indexing
  int matchmethod[n_obs] ; // method indexing
  int matchyears[n_obs]; // year indexing
  int matchsubnat[n_obs]; // subnat indexing 
  vector[n_obs] y; // proportions
  // vector[n_obs] se_prop; // standard errors
  vector[n_years] intercept; // vector 1s
}

parameters {   // The parameters accepted by the model. 
  real alpha_pms[M_count, P_count] ; // expected mean trend
  vector<lower=0>[M_count] sigma_alpha; // variance of mean trend
  vector<lower=0>[M_count] sigma_beta; // world variance
  real beta_c[M_count, C_count]; // overall country mean trend
  vector<lower=0>[P_count] sigmaQ; // sigma RW
  vector[M_count] a_raw[P_count, K]; // variation associated with time
  vector<lower=0>[M_count] tau_a; 
  vector<lower=0>[M_count] sigma_y; 
}

transformed parameters { 
  vector[n_years] z[M_count, P_count]; // latent variable
  matrix<lower=0, upper=1>[S_count, n_years] P[M_count, P_count]; // logit observation
  vector[K] a_k[M_count, P_count]; // variation associated with time

  for(m in 1:M_count){ 
    for(p in 1:P_count){
      a_k[m,p,1] = a_raw[p,1,m];
      for(k in 2:K) {
        a_k[m,p,k] = a_k[m,p,k-1] + a_raw[p,k,m] * tau_a[m]; // penalisiing splines
      }
      z[m,p] = to_vector(alpha_pms[m,p]*intercept) + to_vector(Bik*a_k[m,p]); // Public sector proprtion on logit scale
      P[m, p, 1] = to_row_vector(inv_logit(z[m, p])) ;
      P[m, p, 2] = to_row_vector(1 - P[m, p, 1]) ;
    } // end P loop 
  } // end M loop
}

model { 
  sigma_beta ~ cauchy(0,1);
  sigma_alpha ~ cauchy(0,1); // cross method variance (within a country)
  tau_a ~ cauchy(0,1);
  sigma_y ~ normal(0,2);
  // Priors
  for(m in 1:M_count){ // method loop 1
    for(c in 1:C_count){
      beta_c[m,c] ~ normal(0,sigma_beta[m]);
    }
    for(p in 1:P_count){
      for(k in 1:K) {
        a_raw[p,k,m] ~ normal(0, 1);
      }
    } // end P loop
  }

  
  // Parameter Estimates 
  for(p in 1:P_count) { // country loop
    sigmaQ[p] ~ student_t(2,0,1);
    for(m in 1:M_count){ // method loop 1
      alpha_pms[m,p] ~ normal(beta_c[m, matchcountry[p]],sigma_alpha[m]); // sharing info accross methods within a country so each country public/private sector has an intercept. Tau-alpha is the cross-method variance.
    } // end M loop 1
  } // end P loop
  
  // Likelihood
  for (k in 1:n_obs) {
    y[k] ~ normal(z[matchmethod[k],matchsubnat[k], matchyears[k]],sigma_y[matchmethod[k]]);
  }
}

