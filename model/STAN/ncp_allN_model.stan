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
  }

parameters {   // The parameters accepted by the model. 
  vector[H] delta_k[P_count, M_count]; // variation associated with time
  real<lower=0> sigma_delta; // variance of mean trend
  real<lower=0> sigma_alpha; // variance of mean trend
  real<lower=0> sigma_beta; // variance of mean trend
  vector<lower=0>[M_count] sigma_y; // variance of mean trend
  vector[M_count] alpha_raw[P_count] ; // non-centered parameter for hierarchy
  vector[M_count] beta_c_raw[C_count] ; // expected mean trend
}

transformed parameters { 
  matrix[C_count, M_count] beta_c; // expected mean trend
  matrix[P_count, M_count] alpha_pms; // expected mean trend
  vector[K] beta_k[M_count, P_count]; // spline coefficients
  vector[n_years] z[M_count, P_count]; // latent variable
  matrix<lower=0, upper=1>[S_count, n_years] P[M_count, P_count]; // logit observation
  
  for(m in 1:M_count){ 
    for(c in 1:C_count){
      beta_c[c,m] = sigma_beta*beta_c_raw[c,m];
    }
    for(p in 1:P_count){
      alpha_pms[p, m] = beta_c[matchcountry[p], m] + sigma_alpha*alpha_raw[p,m];
      // Spline coefficients
      beta_k[m,p,kstar[p]] = zero; // set spline coefficient to 0
      for(j in (kstar[p]+1):(K-1)) { // Estimating spline coefficient here
        beta_k[m,p,j] = beta_k[m,p,j-1] + delta_k[p,m, j-1];
      } // after kstar
      for(j in 1:(kstar[p]-1)) { 
        int t = kstar[p] - j;
        beta_k[m,p,t] = beta_k[m,p,t+1] - delta_k[p,m, t];
      } // before kstar
      beta_k[m,p,K] = -(alpha_pms[p, m] + sum(beta_k[m,p,1:(K-1)])); // hard sum to 0 constraint
      // Latent variable
      for(t in 1:n_years) {
        z[m,p,t] = alpha_pms[p,m] +  dot_product(Bik[p, t, 1:K],beta_k[m,p]); // Public sector proprtion on logit scale
      }
      
      // Proportions
      P[m, p, 1] = to_row_vector(inv_logit(z[m, p])) ;
      P[m, p, 2] = to_row_vector(1 - P[m, p, 1]) ;
    } // end P loop 
  } // end M loop
}

model { 
  // Priors
  sigma_delta ~ cauchy(0,1);
  sigma_alpha ~ cauchy(0,1);
  sigma_beta ~ cauchy(0,1);
  sigma_y ~ normal(0,2);

  // Hierarchical estimation of intercept
  for(c in 1:C_count){   // Country intercepts
    beta_c_raw[c] ~ normal(0, 1);
    } // end C loop
    for(p in 1:P_count){
      alpha_raw[p] ~ normal(0,1); // sharing info across methods within a province so each province public/private sector has an intercept.
      for(m in 1:M_count){  
        // sum(beta_k[m,p]) ~ normal(alpha_pms[p,m], inv_sqrt(1 - inv(K)));  // Soft sum-to-0 constraint on spline coefficients 
        delta_k[p,m] ~ normal(0, sigma_delta); // delta are the slopes for logit rates of change in province p, method m, sector s.
      } // end M loop
    } // end P loop

  // Likelihood
  for (k in 1:n_obs) {
    y[k] ~ normal(z[matchmethod[k],matchsubnat[k], matchyears[k]], sigma_y[matchmethod[k]]);
  }
}
