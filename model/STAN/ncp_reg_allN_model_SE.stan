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
  real zero;
  int matchcountry[P_count]; // country indexing
  int matchmethod[n_obs] ; // method indexing
  int matchyears[n_obs]; // year indexing
  int matchsubnat[n_obs]; // subnat indexing 
  vector[n_obs] y; // proportions
  vector[n_obs] se_prop; // standard errors
  real<lower=0> scale_global ; // scale for the half-t prior for tau
  real <lower=0> slab_scale ; # slab scale for the regularized horseshoe
  real <lower=0> slab_df ; # slab degrees of freedom for the regularized horseshoe
  }

parameters {   // The parameters accepted by the model. 
  vector<lower=0>[M_count] sigma_delta; // variance of mean trend
  vector<lower=0>[M_count] sigma_alpha; // variance of mean trend
  vector<lower=0>[M_count] sigma_beta; // variance of mean trend
  vector<lower=0>[M_count] sigma_y; // variance of mean trend
  matrix[M_count, P_count] alpha_raw ; // non-centered parameter for hierarchy
  vector[M_count] beta_c_raw[C_count] ; // expected mean trend
  vector<lower=0>[M_count] tau_delta; 
  vector<lower=0>[M_count] caux;
  vector<lower=0>[H] lambda[P_count, M_count];
  vector[M_count]  logsigma;
}

transformed parameters { 
  vector[H] delta_k[P_count, M_count]; // variation associated with time
  matrix[C_count, M_count] beta_c; // expected mean trend
  matrix[M_count, P_count] alpha_pms; // expected mean trend
  vector[K] beta_k[M_count, P_count]; // spline coefficients
  vector[n_years] z[M_count, P_count]; // latent variable
  matrix[S_count, n_years] P[M_count, P_count]; // logit observation
  vector[H] lstar_sq[P_count, M_count];
  vector<lower=0>[M_count] sigma_tau ;  // noise std
  vector<lower=0>[M_count] cstar; # slab scale
  
  sigma_tau = exp(logsigma);
  
  cstar = slab_scale* sqrt(caux);
  
  for(c in 1:C_count){
    beta_c[c] = to_row_vector(sigma_beta.*beta_c_raw[c]);
  }
  
  for(m in 1:M_count){ 
    for(p in 1:P_count){
      alpha_pms[m,p] = beta_c[matchcountry[p], m] + sigma_alpha[m]*alpha_raw[m,p];
      // Spline coefficients
      beta_k[m,p,kstar[p]] = zero; // set spline coefficient to 0
      lstar_sq[p,m]  = sqrt( cstar[m]^2 * square(lambda[p,m]) ./ (cstar[m]^2 + sigma_tau[m]^2* square(lambda[p,m])) );
      for(j in (kstar[p]+1):K) {
        beta_k[m,p,j] = beta_k[m,p,j-1] + delta_k[p,m, j-1];
      } // after kstar
      for(j in 1:(kstar[p]-1)) { // Estimating spline coefficient here
        int t = kstar[p] - j;
        beta_k[m,p,t] = beta_k[m,p,t+1] - delta_k[p,m, t];
      } // before kstar
      
      delta_k[p,m] = tau_delta[m]*lstar_sq[p,m]; // delta are the slopes for logit rates of change in province p, method m, sector s.

      // Latent variable
      for(t in 1:n_years) {
        z[m,p,t] = alpha_pms[m,p] +  dot_product(Bik[p, t, 1:K],beta_k[m,p]); // Public sector proprtion on logit scale
      }
      
      // Proportions
      P[m, p, 1] = to_row_vector(inv_logit(z[m, p])) ;
      P[m, p, 2] = to_row_vector(1 - P[m, p, 1]) ;
    } // end P loop 
  } // end M loop
}

model { 
  // Priors
  sigma_delta ~ normal(0,2);
  sigma_alpha ~ normal(0,2);
  sigma_beta ~ normal(0,2);
  tau_delta ~ student_t(3 , 0, scale_global*sigma_tau);
  caux ~ inv_gamma (0.5* slab_df, 0.5* slab_df );
  logsigma ~ normal(0, 1);
  
  // Hierarchical estimation of intercept
  for(m in 1:M_count){ 
    for(c in 1:C_count){   // Country intercepts
    beta_c_raw[c,m] ~ normal(0, 1);
    } // end C loop
    for(p in 1:P_count){
      alpha_raw[m,p] ~ normal(0,1); // sharing info across methods within a province so each province public/private sector has an intercept.
      lambda[p,m] ~ cauchy(0,1);
    } // end P loop
  } // end M loop

  // Likelihood
  for (k in 1:n_obs) {
    y[k] ~ normal(z[matchmethod[k],matchsubnat[k], matchyears[k]], se_prop[k]);
  }
}
