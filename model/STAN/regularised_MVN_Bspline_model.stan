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

}

parameters {
  vector[num_knots] a_raw[P_count, M_count];
  vector[num_knots] lambda_araw[P_count, M_count]; // allows movement away from 0
  vector[M_count] a0[P_count];  // intercept
  vector<lower=0>[M_count] sigma_y; // data variance
  vector<lower=0>[M_count] sigma_a0; // intercept variance
  vector<lower=0>[M_count] tau; // spline variance
  vector<lower=0>[M_count] tau_araw; // regulaisation parameter
  cholesky_factor_corr[M_count] sigmaa0_Omega; // prior correlation
  vector<lower=0>[M_count] sigmaa0_tau;  // prior scale
  vector[M_count] beta_c;
  vector[M_count] logsigma;
  real c_sq;
}

transformed parameters {
  vector[num_knots] a[P_count, M_count];
  vector[num_years] Y_hat[P_count, M_count];
  matrix<lower=0, upper=1>[S_count, num_years] P[P_count, M_count]; // logit observation
  vector[num_knots] lstar_sq[P_count, M_count];
  vector<lower=0>[M_count] sigma ;  // noise std
  cholesky_factor_cov[M_count] L_Sigma_a0; // cholesky variance of province-level mean trend

  sigma = exp(logsigma);
  
  L_Sigma_a0 = diag_pre_multiply(sigmaa0_tau, sigmaa0_Omega);

  for(p in 1:P_count) {
    for(m in 1:M_count) {
      for(k in 1:num_knots){
        lstar_sq[p,m,k] = (c_sq*pow(lambda_araw[p,m,k],2))/(c_sq+pow(tau_araw[m],2)*pow(lambda_araw[p,m,k],2));
      }
      
      a[p,m,1] = a_raw[p,m,1];
      for (i in 2:num_knots) {
        a[p,m,i] = a[p,m, i-1] + a_raw[p,m,i]*tau[m]; // penalising splines
      }
      Y_hat[p,m,1:num_years] = a0[p,m]*to_vector(X) + to_vector(B'*a[p,m, 1:num_knots]);
      P[p, m, 1, 1:num_years] = to_row_vector(inv_logit(Y_hat[p,m,1:num_years])) ;
      P[p, m, 2, 1:num_years] = to_row_vector(1 - P[p, m, 1, 1:num_years]) ;
    }
  }
}

model {
  // Priors
  tau ~ normal(0, 2);
  sigma_y ~ normal(0, 2);
  sigma_a0 ~ normal(0, 1);
  beta_c ~ normal(0, 2);
  c_sq ~ inv_gamma(2, 8);
  sigmaa0_tau ~ normal(0,2);
  sigmaa0_Omega ~ lkj_corr_cholesky(1);
  
  for(m in 1:M_count){
    tau_araw[m] ~ student_t(3 , 0, scale_global*sigma[m]);
  }
  
  for(p in 1:P_count){
    a0[p] ~ multi_normal_cholesky(beta_c, L_Sigma_a0);
    for(m in 1:M_count){
      for(k in 1:num_knots) {
        lambda_araw[p,m,k] ~ cauchy(0,1);
        a_raw[p,m,k] ~ normal(0, tau_araw[m]*sqrt(lstar_sq[p,m,k])); 
      }
    }
  }
  
  // Likelihood
  for (k in 1:num_data) {
    Y[k] ~ normal(Y_hat[matchsubnat[k], matchmethod[k], matchyears[k]], sigma_y[matchmethod[k]]);
  }
}

