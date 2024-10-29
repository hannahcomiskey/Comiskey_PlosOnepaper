data {  
  int<lower=1> n_years; // Number of years
  int<lower=1> n_obs; // Number of observations
  int<lower=1> H; // Number of knots
  int<lower=1> K; // Number of spline cofficients (H+1)
  int<lower=1> P_count; // Number of provinces
  int<lower=1> C_count; // Number of countries
  int<lower=1> M_count; // Number of methods
  int<lower=1> S_count; // Number of sectors
  array[P_count] int<lower=1, upper=K> kstar; // Spline index K star for estimation
  vector[K] Bik[P_count, n_years]; // Basis functions
  int zero;
  int matchcountry[P_count]; // country indexing
  int matchmethod[n_obs] ; // method indexing
  int matchyears[n_obs]; // year indexing
  int matchsubnat[n_obs]; // subnat indexing 
  vector[n_obs] y; // proportions
}

parameters {   // The parameters accepted by the model. 
  vector<lower=-9, upper=9>[M_count] alpha_pms[P_count]; // expected mean trend
  cholesky_factor_corr[M_count] sigmaalpha_Omega; // prior correlation
  vector<lower=0>[M_count] sigmaalpha_tau;  // prior scale
  vector<lower=-9, upper=9>[M_count] beta_c[C_count]; // overall country mean trend
  cholesky_factor_corr[M_count] sigmabeta_Omega; // prior correlation
  vector<lower=0>[M_count] sigmabeta_tau;  // prior scale
  vector[H] delta_k[P_count, M_count]; // variation associated with time
  vector<lower=0>[M_count] sigma_delta; // variance of mean trend
}

transformed parameters { 
  vector[K] beta_k[M_count, P_count]; // spline coefficients
  vector[n_years] z[M_count, P_count]; // latent variable
  matrix[S_count, n_years] P[M_count, P_count]; // logit observation
  cholesky_factor_cov[M_count] L_Sigma_beta; // cholesky variance of country-level mean trend
  cholesky_factor_cov[M_count] L_Sigma_alpha; // cholesky variance of province-level mean trend

  L_Sigma_beta = diag_pre_multiply(sigmabeta_tau, sigmabeta_Omega);
  L_Sigma_alpha = diag_pre_multiply(sigmaalpha_tau, sigmaalpha_Omega);
  
  for(m in 1:M_count){ 
    for(p in 1:P_count){
      
      // Spline coefficients
      beta_k[m,p,kstar[p]] = zero; // set spline coefficient to 0
      for(j in (kstar[p]+1):K) {
        beta_k[m,p,j] = beta_k[m,p,j-1] + delta_k[p, m, j-1];
      } // after kstar
      for(j in 1:(kstar[p]-1)) { // Estimating spline coefficient here
        int t = kstar[p] - j;
        beta_k[m,p,t] = beta_k[m,p,t+1] - delta_k[p, m, t];
      } // before kstar
      
      // Latent variable
      for(t in 1:n_years) {
        z[m,p,t] = alpha_pms[p,m] + dot_product(Bik[p,t],beta_k[m,p]); // Public sector proprtion on logit scale
      }
      
      // Proportions
      P[m, p, 1] = to_row_vector(inv_logit(z[m, p])) ;
      P[m, p, 2] = to_row_vector(1 - P[m, p, 1]) ;
    } // end P loop 
  } // end M loop
}

model { 
  // Priors
  sigmabeta_tau ~ normal(0,2);
  sigmabeta_Omega ~ lkj_corr_cholesky(1);
  sigmaalpha_tau ~ normal(0,2);
  sigmaalpha_Omega ~ lkj_corr_cholesky(1);
  sigma_delta ~ normal(0,2); // cross-country variance (within a method)
 
  // Hierarchical estimation of intercept
  for(c in 1:C_count){   // Country intercepts
      beta_c[c] ~ multi_normal_cholesky(rep_vector(0, M_count), L_Sigma_beta);
    } // end C loop
  for(p in 1:P_count){
    alpha_pms[p] ~ multi_normal_cholesky(beta_c[matchcountry[p]], L_Sigma_alpha); // sharing info across methods within a province so each province public/private sector has an intercept.
    for(m in 1:M_count){
      for(h in 1:H){
        delta_k[p,m,h] ~ normal(0, sigma_delta[m]); // delta are the slopes for logit rates of change in province p, method m, sector s.
      } // end H loop
    } // end M loop
  } // end P loop

  // Likelihood
  for (k in 1:n_obs) {
    y[k] ~ normal(z[matchmethod[k],matchsubnat[k], matchyears[k]], 1);
  }
}


