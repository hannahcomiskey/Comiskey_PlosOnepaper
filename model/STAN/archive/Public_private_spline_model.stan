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
  int kstar[P_count]; // indexing K star for estimation
  vector[P_count] splinestart; // vector of 0s
  int matchcountry[P_count]; // country indexing
  int matchmethod[n_obs] ; // method indexing
  int matchyears[n_obs]; // year indexing
  int matchsubnat[n_obs]; // subnat indexing 
  vector[n_obs] y; // proportions
  vector[n_obs] se_prop; // standard errors
  vector[n_years] intercept; // vector 1s
  vector[M_count] sigmadelta_mu; // vector of 0s
}

parameters {   // The parameters accepted by the model. 
  real alpha_pms[M_count, P_count] ; // expected mean trend
  vector<lower=0>[C_count] sigma_alpha; // variance of mean trend
  real beta_world; // world mean trend for sectors
  real <lower=0> sigma_beta; // world variance
  vector[C_count] beta_c; // overall country mean trend
  vector[M_count] delta_k[C_count, H]; // variation associated with time 
  corr_matrix[M_count] sigmadelta_Omega; // prior correlation
  vector<lower=0>[M_count] sigmadelta_tau;  // prior scale
}

transformed parameters { 
  vector[K] beta_k[M_count, P_count]; // spline coefficients
  vector[n_years] z[M_count, P_count]; // latent variable
  for(m in 1:M_count){ 
    for(p in 1:P_count){
      beta_k[m,p,kstar[p]] = splinestart[p];
      for(j in 1:(kstar[p]-1)) { // Estimating spline coefficient here
         beta_k[m,p,(kstar[p] - j)] = beta_k[m,p,(kstar[p] - j)+1] - delta_k[p,(kstar[p] - j),m];
      } // end K1 loop (before kstar)
      for(j in (kstar[p]+1):K) {
        beta_k[m,p,j] = beta_k[m,p,(j-1)] + delta_k[p,(j-1),m];
      } // end K2 loop (after kstar)
      z[m,p] = to_vector(alpha_pms[m,p]*intercept) + to_vector(Bik[p]*beta_k[m,p]); // Public sector proprtion on logit scale
    } // end P loop 
  } // end M loop
}

model { 
  // Priors
  sigmadelta_tau ~ cauchy(0, 2.5);
  sigmadelta_Omega ~ lkj_corr(1);
  sigma_beta ~ student_t(3,0,1);   // cross country variance
  beta_world ~ normal(0,10);  // world intercept
  for(c in 1:C_count){   // Country intercepts 
    beta_c[c]~ normal(beta_world,sigma_beta);
    sigma_alpha[c] ~ student_t(1,0,1); // cross method variance (within a country)
  }
  
  // Parameter Estimates 
  for(p in 1:P_count) { // country loop
    for(m in 1:M_count){ // method loop 1
      alpha_pms[m,p] ~ normal(beta_c[matchcountry[p]],sigma_alpha[matchcountry[p]]); // sharing info accross methods within a country so each country public/private sector has an intercept. Tau-alpha is the cross-method variance.
    } // end M loop 1
    for(j in 1:H){ 
      delta_k[p,j] ~ multi_normal(sigmadelta_mu, quad_form_diag(sigmadelta_Omega, sigmadelta_tau)); // delta are the slopes for logit rates of change in province p, method m, sector s.
    } // end H loop
  } // end P loop
  
  // Likelihood
  for (k in 1:n_obs) {
    y[k] ~ normal(z[matchmethod[k],matchsubnat[k], matchyears[k]], se_prop[k]);
  }
}

