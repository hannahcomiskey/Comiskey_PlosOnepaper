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
  real<lower=0> scale_global ; // scale for the half-t prior for tau
}

parameters {   // The parameters accepted by the model. 
  // corr_matrix[M_count] sigmaalpha_Omega; // prior correlation
  // vector<lower=0>[M_count] sigmaalpha_tau;  // prior scale
  vector<lower=0>[M_count] sigma_alpha;
  corr_matrix[M_count] sigmabeta_Omega; // prior correlation
  vector<lower=0>[M_count] sigmabeta_tau;  // prior scale
  vector[H] delta_k[P_count, M_count]; // variation associated with time
  vector<lower=0>[M_count] sigma_y; // data variance
  vector[M_count] delta_betac[C_count]; // variation for beta parameters
  vector[M_count] alpha_raw[P_count]; // variation for beta parameters
  real<lower=0> tau; 
  real c_sq;
  vector<lower=0>[H] lambda[P_count, M_count];
  real logsigma;

}

transformed parameters { 
  vector[K] beta_k[M_count, P_count]; // spline coefficients
  vector[n_years] z[M_count, P_count]; // latent variable
  matrix[S_count, n_years] P[M_count, P_count]; // logit observation
  vector[M_count] beta_c[C_count]; // overall country mean trend
  vector[M_count] alpha_pms[P_count]; // expected mean trend
  vector[H] lstar_sq[P_count, M_count];
  real <lower=0> sigma ;  // noise std
  
  sigma = exp(logsigma);
  
  for(c in 1:C_count){
    for(m in 1:M_count){ 
    beta_c[c, 1:M_count] =  quad_form_diag(sigmabeta_Omega, sigmabeta_tau)*delta_betac[c, 1:M_count];
    // beta_c[c, m] =   sigma_beta[m]*beta_raw[c,m];
    }
  }
  for(p in 1:P_count){
    // alpha_pms[p,1:M_count] = beta_c[matchcountry[p], 1:M_count] +  quad_form_diag(sigmaalpha_Omega, sigmaalpha_tau)*delta_alpha[p, 1:M_count];  
    for(m in 1:M_count){ 
      alpha_pms[p,m] = beta_c[matchcountry[p], m] + sigma_alpha[m]*alpha_raw[p,m];
      for(h in 1:H){
        lstar_sq[p,m,h] = (c_sq*pow(lambda[p,m,h],2))/(c_sq+pow(tau,2)*pow(lambda[p,m,h],2));
      }
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
  sigma_alpha ~ normal(0,2);
  sigmabeta_tau ~ normal(0,1);
  sigmabeta_Omega ~ lkj_corr(2);
  // sigmaalpha_tau ~ normal(0,1);
  // sigmaalpha_Omega ~ lkj_corr(2);
  sigma_y ~ normal(0, 2);
  c_sq ~ inv_gamma(2, 8);
  tau ~ student_t(3 , 0, scale_global*sigma);
  
  // Hierarchical estimation of intercept
  for(c in 1:C_count){   // Country intercepts
    delta_betac[c] ~ normal(0,1); 
  } // end C loop
  for(p in 1:P_count){
    alpha_raw[p] ~ normal(0,1); // sharing info across methods within a province so each province public/private sector has an intercept.
    for(m in 1:M_count){
      for(h in 1:H){
        lambda[p,m,h] ~ cauchy(0,1);
        delta_k[p,m,h] ~ normal(0, tau*sqrt(lstar_sq[p,m,h])); // delta are the slopes for logit rates of change in province p, method m, sector s.
      } // end H loop
    } // end M loop
  } // end P loop

  // Likelihood
  for (k in 1:n_obs) {
    y[k] ~ normal(z[matchmethod[k],matchsubnat[k], matchyears[k]], sigma_y[matchmethod[k]]);
  }
}


