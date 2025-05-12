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
  real<lower=1> nu_global ; // degrees of freedom for the half-t prior for tau
  real<lower=1> nu_local ; // degrees of freedom for the half-t priors for lambdas
  real<lower=0> scale_global ; // scale for the half-t prior for tau
  real<lower=0> slab_scale ; # slab scale for the regularized horseshoe
  real<lower=0> slab_df ; # slab degrees of freedom for the regularized horseshoe
  }

parameters {   // The parameters accepted by the model.
  vector[H] delta_k[P_count, M_count]; // variation associated with time
  vector<lower=0>[M_count] sigma_delta; // variance of mean trend
  vector<lower=0>[M_count] sigma_alpha; // variance of mean trend
  vector<lower=0>[M_count] sigma_beta; // variance of mean trend
  vector<lower=0>[M_count] sigma_y;
  vector[M_count] alpha_raw[P_count]; // non-centered parameter for hierarchy
  vector[M_count] beta_c_raw[C_count] ; // expected mean trend
  real<lower=0> caux;
  real logsigma;
  real<lower=0> aux1_global;
  real<lower=0> aux2_global;
  vector<lower=0>[H] aux1_local[P_count, M_count];
  vector<lower=0>[H] aux2_local[P_count, M_count];
}

transformed parameters { 
  vector[M_count] beta_c[C_count]; // expected mean trend
  vector[M_count] alpha_pms[P_count]; // expected mean trend
  vector[K] beta_k[M_count, P_count]; // spline coefficients
  vector[n_years] z[M_count, P_count]; // latent variable
  matrix[S_count, n_years] P[M_count, P_count]; // logit observation
  vector<lower=0>[H] lambda_tilde[P_count, M_count]; // 'truncated' local shrinkage parameter
  vector<lower=0>[H] lambda[P_count, M_count]; // local shrinkage parameter
  real<lower=0> sigma_tau;  // noise std
  real<lower=0> cstar; # slab scale
  real<lower=0> tau_delta; // global shrinkage parameter

  sigma_tau = exp(logsigma);
  tau_delta = aux1_global*sqrt(aux2_global)*scale_global*sigma_tau;
  cstar = slab_scale*sqrt(caux);
  
  for(c in 1:C_count){
    beta_c[c] = sigma_beta.*beta_c_raw[c];
  }
  
  for(m in 1:M_count){ 
    for(p in 1:P_count){
      alpha_pms[p,m] = beta_c[matchcountry[p], m] + sigma_alpha[m]*alpha_raw[p,m];
      // Spline coefficients
      beta_k[m,p,kstar[p]] = zero; // set spline coefficient to 0
      for(j in (kstar[p]+1):K) {
        beta_k[m,p,j] = beta_k[m,p,j-1] + delta_k[p,m, j-1];
      } // after kstar
      for(j in 1:(kstar[p]-1)) { // Estimating spline coefficient here
        int t = kstar[p] - j;
        beta_k[m,p,t] = beta_k[m,p,t+1] - delta_k[p,m, t];
      } // before kstar
      lambda[p,m] = aux1_local[p,m].*sqrt(aux2_local[p,m]);
      lambda_tilde[p,m]  = sqrt((pow(cstar,2)*pow(lambda[p,m],2)) ./ (pow(cstar,2) + pow(sigma_tau,2)*pow(lambda[p,m], 2)));
      // delta_k[p,m,h] = tau_delta[m]*lambda_tilde[p,m,h]; // delta are the slopes for logit rates of change in province p, method m, sector s.

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
  sigma_delta ~ normal(0, 2);
  sigma_alpha ~ normal(0, 2);
  sigma_beta ~ normal(0, 2);
  caux ~ inv_gamma(0.5*slab_df, 0.5*slab_df);
  logsigma ~ normal(0, 2);
  aux1_global ~ normal(0, 1);
  aux2_global ~ inv_gamma(0.5*nu_global, 0.5*nu_global);
  sigma_y ~ normal(0, 2);
  
  // Hierarchical estimation of intercept
  for(c in 1:C_count){   // Country intercepts
    beta_c_raw[c] ~ normal(0, 1);
  } // end C loop
  
  for(p in 1:P_count){
    alpha_raw[p] ~ normal(0, 1); // sharing info across methods within a province so each province public/private sector has an intercept.
  } // end P loop
  
  for(m in 1:M_count) {
    // tau_delta[m] ~ student_t(3, 0, scale_global*sigma_tau[m]);
    for(p in 1:P_count){  
      // lambda[p,m] ~ cauchy(0,1);
      sum(beta_k[m,p]) ~ normal(alpha_pms[p,m], inv_sqrt(1 - inv(K)));  // Sum-to-0 constraint on spline coefficients 
      for(h in 1:H) {
        delta_k[p,m,h] ~ normal(0, tau_delta*lambda_tilde[p,m,h]);
      }
      aux1_local[p,m] ~ normal(0, 1);
      aux2_local[p,m] ~ inv_gamma(0.5*nu_local, 0.5*nu_local);
    }
  }


  // Likelihood
  for (k in 1:n_obs) {
    y[k] ~ normal(z[matchmethod[k],matchsubnat[k], matchyears[k]], sigma_y[matchmethod[k]]);
  }
}

