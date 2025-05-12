data {  
  int<lower=1> n_years; // Number of years
  int<lower=1> n_obs; // Number of observations
  int<lower=1> H; // Number of knots
  int<lower=1> K; // Number of spline cofficients (H+1)
  int<lower=1> P_count; // Number of provinces
  int<lower=1> C_count; // Number of countries
  int<lower=1> M_count; // Number of methods
  int<lower=1> S_count; // Number of sectors
  matrix[n_years, H] Zih; // Basis functions
  vector[n_years] intercept; // vector 1s
  vector[M_count] beta_mu; // vector 0s
  vector[M_count] delta_mu; // vector 0s
  int matchcountry[P_count]; // country indexing
  int matchmethod[n_obs] ; // method indexing
  int matchyears[n_obs]; // year indexing
  int matchsubnat[n_obs]; // subnat indexing 
  vector[n_obs] y; // proportions
  vector<lower=0>[n_obs] se_prop; // standard errors
}

parameters {   // The parameters accepted by the model. 
  real alpha_pms[M_count, P_count] ; // expected mean trend
  vector<lower=0>[M_count] sigma_alpha; // variance of mean trend
  vector<lower=0>[M_count] sigma_delta; // variance of mean trend
  vector[M_count] beta_c[C_count]; // overall country mean trend
  cholesky_factor_corr[M_count] sigmabeta_Omega; // prior correlation
  vector<lower=0>[M_count] sigmabeta_tau;  // prior scale
  vector[H-1] delta_k_raw[M_count, P_count]; // variation associated with time
}

transformed parameters { 
  vector[H] delta_k[M_count, P_count]; // variation associated with time
  vector[n_years] z[M_count, P_count]; // latent variable
  matrix<lower=0, upper=1>[S_count, n_years] P[M_count, P_count]; // logit observation
  matrix[M_count, M_count] L_Sigma; // cholesky decomp. of covariance

  L_Sigma = diag_pre_multiply(sigmabeta_tau, sigmabeta_Omega);

  for(m in 1:M_count){ 
    for(p in 1:P_count){
      delta_k[m,p] = append_row(delta_k_raw[m,p], -sum(delta_k_raw[m,p])); // sum 0 constraint
      z[m,p] = to_vector(alpha_pms[m,p]*intercept) + to_vector(Zih*delta_k[m,p]); // Public sector proprtion on logit scale
      P[m, p, 1] = to_row_vector(inv_logit(z[m, p])) ;
      P[m, p, 2] = to_row_vector(1 - P[m, p, 1]) ;
    } // end P loop 
  } // end M loop
}

model { 
  sigmabeta_tau ~ cauchy(0, 2.5);
  sigmabeta_Omega ~ lkj_corr_cholesky(4);

  // Priors
  for(m in 1:M_count){ 
    sigma_alpha[m] ~ cauchy(0, 2.5); // cross-country variance (within a method)
    sigma_delta[m] ~ cauchy(0, 2.5); // cross-country variance (within a method)
    for(p in 1:P_count){ 
      alpha_pms[m,p] ~ normal(beta_c[matchcountry[p], m],sigma_alpha[m]); // sharing info accross methods within a country so each country public/private sector has an intercept. Tau-alpha is the cross-method variance.
      for(h in 1:(H-1)){
        delta_k_raw[m,p,h] ~ normal(0, sigma_delta[m]); // multi_normal(delta_mu, quad_form_diag(sigmadelta_Omega, sigmadelta_tau));  delta are the slopes for logit rates of change in province p, method m, sector s.
      } // end H loop
    } // end P loop
  } // end M loop
  for(c in 1:C_count){   // Country intercepts
    beta_c[c] ~ multi_normal_cholesky(beta_mu, L_Sigma);
  } // end C loop
 
  // Likelihood
  for (k in 1:n_obs) {
    y[k] ~ normal(z[matchmethod[k],matchsubnat[k], matchyears[k]], se_prop[k]);
  }
}

generated quantities {
  vector[n_obs] y_tilde;
  for (k in 1:n_obs) {
    y_tilde[k] = normal_rng(z[matchmethod[k],matchsubnat[k], matchyears[k]], se_prop[k]);
  }
}

