data {  
  int<lower=0> n_years; // Number of years
  int<lower=0> n_obs; // Number of observations
  int<lower=0> H; // Number of knots
  int<lower=0> K; // Number of spline cofficients (H+1)
  int<lower=0> P_count; // Number of provinces
  int<lower=0> C_count; // Number of countries
  int<lower=0> M_count; // Number of methods
  int<lower=0> S_count; // Number of sectors
  matrix[n_years, K] Bik[P_count]; // Basis functions
  vector[n_years] intercept; // vector 1s
  vector[M_count] beta_mu; // vector 0s
  int matchcountry[P_count]; // country indexing
  int matchmethod[n_obs] ; // method indexing
  int matchyears[n_obs]; // year indexing
  int matchsubnat[n_obs]; // subnat indexing 
  vector[n_obs] y; // proportions
  vector[n_obs] se_prop; // standard errors
  // matrix[n_obs, S_count] Y; // proportions
  // matrix[S_count, S_count] Sigma_Y[n_obs]; // variance-covariance matrix
}

parameters {   // The parameters accepted by the model. 
  real alpha_pms[M_count, P_count] ; // expected mean trend
  vector<lower=0>[C_count] sigma_alpha; // variance of mean trend
  matrix[C_count, M_count] beta_c; // overall country mean trend
  vector[K] beta_k[M_count, P_count]; // spline coefficients
  corr_matrix[M_count] sigmabeta_Omega; // prior correlation
  vector<lower=0>[M_count] sigmabeta_tau;  // prior scale
}

transformed parameters { 
  vector[n_years] z[M_count, P_count]; // latent variable
  matrix[S_count, n_years] logit_estimates[M_count, P_count]; // logit observations
  matrix[S_count, n_years] P[M_count, P_count]; // logit observations
  for(m in 1:M_count){ 
    for(p in 1:P_count){
      z[m,p] = to_vector(alpha_pms[m,p]*intercept) + to_vector(Bik[p]*beta_k[m,p]); // Public sector proprtion on logit scale
      P[m, p, 1] = to_row_vector(inv_logit(z[m, p])) ;
      P[m, p, 2] = to_row_vector(1 - P[m, p, 1]) ;
      logit_estimates[m, p, 1] = to_row_vector(z[m,p]) ; // logit public observations
      logit_estimates[m, p, 2] = to_row_vector(logit(P[m, p, 2])) ; // logit private observations
    } // end P loop 
  } // end M loop
}

model { 
  // Priors
  sigmabeta_tau ~ cauchy(0, 2.5);
  sigmabeta_Omega ~ lkj_corr(1);
  
  // Parameter Estimates 
  for(c in 1:C_count){   // Country intercepts 
    beta_c[c] ~ multi_normal(beta_mu, quad_form_diag(sigmabeta_Omega, sigmabeta_tau));
  } // end C loop
  for(m in 1:M_count){ 
    sigma_alpha[m] ~ normal(0,1); // cross-country variance (within a method)
    for(p in 1:P_count){ 
      alpha_pms[m,p] ~ normal(beta_c[matchcountry[p], m],sigma_alpha[m]); // sharing info accross methods within a country so each country public/private sector has an intercept. Tau-alpha is the cross-method variance.
      for(k in 1:(K-1)){
        beta_k[m,p,k] ~ normal(0, 1);
      }
      sum(beta_k[m,p]) ~ normal(0, 0.001*S_count); // soft sum 0 constraint
    } // end P loop 
  } // end M loop

  // Likelihood
  for (k in 1:n_obs) {
    y[k] ~ normal(z[matchmethod[k],matchsubnat[k], matchyears[k]], se_prop[k]);
    // Y[k, 1:2] ~ multi_normal(to_vector(logit_estimates[matchmethod[k], matchsubnat[k]]), Sigma_Y[k]);
  }
}


