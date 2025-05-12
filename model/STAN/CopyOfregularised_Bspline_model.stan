data {
  int num_data;             // number of data points
  int num_knots;            // num of knots
  int num_years;
  vector[num_data] Y;
  vector[num_years] X;
  matrix[num_knots, num_years] B;
  int M_count;                    // Number of methods
  int P_count;                    // Number of provinces
  int S_count;
  int matchmethod[num_data] ; // method indexing
  int matchyears[num_data]; // year indexing
  int matchsubnat[num_data]; // subnat indexing 
  real<lower=0> scale_global ; // scale for the half -t prior for tau
  real<lower=0> slab_scale ; # slab scale for the regularized horseshoe
  real<lower=0> slab_df ; # slab degrees of freedom for the regularized horseshoe
}

parameters {
  vector<lower=0>[num_knots] lambda[P_count, M_count]; // local shrinkage parameter
  vector[M_count] a0_raw[P_count];  // intercept
  vector<lower=0>[M_count] sigma_y; // data variance
  vector<lower=0>[M_count] sigma_a0; // intercept variance
  vector<lower=0>[M_count] sigma_tau;  // noise std
  vector<lower=0>[M_count] tau_a; // spline variance
  vector[M_count] beta_c;
  vector<lower=0>[M_count] caux;
  vector<lower=0>[M_count] tau_a0; // global shrinkage parameter
}

transformed parameters {
  vector[num_knots] a[P_count, M_count];
  vector[num_knots] a_raw[P_count, M_count];
  vector[M_count] a0[P_count];  // intercept
  vector[num_years] Y_hat[P_count, M_count];
  matrix<lower=0, upper=1>[S_count, num_years] P[P_count, M_count]; // logit observation
  vector[num_knots] lambda_tilde[P_count, M_count];
  vector<lower=0>[M_count] cstar; # slab scale
  
  for(m in 1:M_count) {
    cstar[m] = slab_scale*sqrt(caux[m]);
    for(p in 1:P_count) {
      a0[p,m] = beta_c[m] + sigma_a0[m]*a0_raw[p,m]; // sharing info across methods within a province so each province public/private sector has an intercept.
      for(k in 1:num_knots){
        lambda_tilde[p,m,k]  = sqrt((pow(cstar[m],2)*pow(lambda[p,m,k],2))/ (pow(cstar[m],2) + pow(sigma_tau[m],2)*pow(lambda[p,m,k], 2)));
        a_raw[p,m,k] = tau_a0[m]*lambda_tilde[p,m,k];
      }
      a[p,m,1] = a_raw[p,m,1];
      for (i in 2:num_knots) {
        a[p,m,i] = a[p,m, i-1] + a_raw[p,m,i]*tau_a[m]; // penalising splines
      }
      Y_hat[p,m,1:num_years] = a0[p,m]*to_vector(X) + to_vector(B'*a[p,m, 1:num_knots]);
      P[p, m, 1, 1:num_years] = to_row_vector(inv_logit(Y_hat[p,m,1:num_years])) ;
      P[p, m, 2, 1:num_years] = to_row_vector(1 - P[p, m, 1, 1:num_years]) ;
    }
  }
}

model {
  // Priors
  tau_a ~ normal(0, 2);
  sigma_y ~ normal(0, 2);
  sigma_a0 ~ normal(0, 1);
  beta_c ~ normal(0, 2);
  caux ~ inv_gamma(0.5*slab_df, 0.5*slab_df);
  sigma_tau ~ cauchy(0,1);
  
  for(m in 1:M_count){
    tau_a0[m] ~ student_t(3, 0, scale_global*sigma_tau[m]);
    for(p in 1:P_count){  
      lambda[p,m] ~ cauchy(0,1);
      a0_raw[p,m] ~ normal(0,1);
      sum(a[p,m]) ~ normal(0, inv_sqrt(1 - inv(num_knots))); // soft sum-to-0 constraint
    }
  }

  // Likelihood
  for (k in 1:num_data) {
    Y[k] ~ normal(Y_hat[matchsubnat[k], matchmethod[k], matchyears[k]], sigma_y[matchmethod[k]]);
  }
}

