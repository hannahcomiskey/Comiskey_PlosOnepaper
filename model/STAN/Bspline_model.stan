data {
  int num_data;             // number of data points
  int num_knots;            // num of knots
  int num_years;
  vector[num_data] Y;
  vector[num_years] X;
  matrix[num_years, num_knots] B;
  int M_count;                    // Number of methods
  int P_count;                    // Number of provinces
  int S_count;
  int matchmethod[num_data] ; // method indexing
  int matchyears[num_data]; // year indexing
  int matchsubnat[num_data]; // subnat indexing 
}

parameters {
  vector[num_knots] a_raw[P_count, M_count];
  vector<lower=0>[M_count] sigma_y; // data variance
  vector<lower=0>[M_count] sigma_a0; // intercept variance
  vector<lower=0>[M_count] tau; // spline variance
  matrix[P_count, M_count] a0_raw;
  vector[M_count] beta_c;
}

transformed parameters {
  vector[num_knots] a[P_count, M_count];
  vector[num_years] Y_hat[P_count, M_count];
  matrix<lower=0, upper=1>[S_count, num_years] P[P_count, M_count]; // logit observation
  matrix[P_count, M_count] a0;  // intercept
  
  for(p in 1:P_count) {
    for(m in 1:M_count) {
      a0[p,m] = beta_c[m] + sigma_a0[m]*a0_raw[p,m]; // NCP parameterization
      a[p,m,1] = a_raw[p,m,1];
      for (i in 2:(num_knots-1)) {
        a[p,m,i] = a[p,m, i-1] + a_raw[p,m,i]*tau[m]; // penalising splines
      }
      a[p,m,num_knots] = -sum(a[p,m,1:(num_knots-1)]); # hard sum to 0 constraint
      Y_hat[p,m,1:num_years] = a0[p,m]*to_vector(X) + to_vector(B*a[p,m, 1:num_knots]);
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
  
  for(p in 1:P_count){
    //beta_p ~ normal(0, 2);
    for(m in 1:M_count){
      a0_raw[p,m] ~ normal(0, 1); 
      a_raw[p,m] ~ normal(0, 1);
    }
  }

  // Likelihood
  for (k in 1:num_data) {
    Y[k] ~ normal(Y_hat[matchsubnat[k], matchmethod[k], matchyears[k]], sigma_y[matchmethod[k]]);
  }
}

