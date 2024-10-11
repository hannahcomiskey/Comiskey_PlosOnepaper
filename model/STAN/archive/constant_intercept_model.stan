data{  // The input data
  int<lower=0> n_years; // Number of years
  int<lower=0> n_obs; // Number of observations
  int<lower=0> H; // Number of knots
  int<lower=0> K; // Number of spline cofficients (H+1)
  int<lower=0> P_count; // Number of provinces
  int<lower=0> C_count; // Number of countries
  int<lower=0> M_count; // Number of methods
  int<lower = 1> kstar[P_count]; // Spline index K star for estimation
  int matchcountry[P_count]; // country indexing
  int matchmethod[n_obs] ; // method indexing
  int matchyears[n_obs]; // year indexing
  int matchsubnat[n_obs]; // subnat indexing 
  vector[n_obs] y; // proportions
  vector[n_obs] se_prop; // standard errors
}

parameters {   // The parameters accepted by the model. 
  real alpha_pms[M_count, P_count] ; // expected mean trend
  real x[P_count, n_years] ; // Random walk
  vector<lower=0>[C_count] sigma_alpha; // variance of mean trend
  real beta_world; // world mean trend for sectors
  real <lower=0> sigma_beta; // world variance
  vector[C_count] beta_c; // overall country mean trend
  vector<lower=0>[P_count] sigmaQ; // sigma RW
}

transformed parameters { 
  real z[M_count, P_count, n_years]; // latent variable
  for(m in 1:M_count){ 
    for(p in 1:P_count){
      for (t in 1:n_years) { // logit-public
        z[m,p,t] = alpha_pms[m,p] + x[p,t]; // Public sector proprtion on logit scale
      } // end t loop
    } // end P loop 
  } // end M loop
}

model { 
  // Priors
  sigma_beta ~ student_t(3,0,1);   // cross country variance
  beta_world ~ normal(0,10);  // world intercept
  for(c in 1:C_count){   // Country intercepts 
    beta_c[c]~ normal(beta_world,sigma_beta);
    sigma_alpha[c] ~ student_t(1,0,1); // cross method variance (within a country)
  }
  
  // Parameter Estimates 
  for(p in 1:P_count) { // country loop
    sigmaQ[p] ~ student_t(2,0,1);
    x[p,1] ~ normal(0, 1);
    for (t in 2:n_years) {
      x[p,t] ~ normal(x[p,t-1], sigmaQ[p]);
    }
    for(m in 1:M_count){ // method loop 1
      alpha_pms[m,p] ~ normal(beta_c[matchcountry[p]],sigma_alpha[matchcountry[p]]); // sharing info accross methods within a country so each country public/private sector has an intercept. Tau-alpha is the cross-method variance.
    } // end M loop 1
  } // end P loop
  
  // Likelihood
  for (k in 1:n_obs) {
    y[k] ~ normal(z[matchmethod[k],matchsubnat[k], matchyears[k]], se_prop[k]);
  }
}

