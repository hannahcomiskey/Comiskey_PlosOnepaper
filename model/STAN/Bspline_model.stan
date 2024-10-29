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
}

parameters {
  vector[num_knots] a_raw[P_count, M_count];
  vector[M_count] a0;  // intercept
  vector<lower=0>[M_count] sigma_y; // data variance
  vector<lower=0>[M_count] tau; // spline variance
}

transformed parameters {
  vector[num_knots] a[P_count, M_count];
  vector[num_years] Y_hat[P_count, M_count];
  matrix<lower=0, upper=1>[S_count, num_years] P[P_count, M_count]; // logit observation

  for(p in 1:P_count) {
    for(m in 1:M_count) {
      a[p,m,1] = a_raw[p,m,1];
      for (i in 2:num_knots) {
        a[p,m,i] = a[p,m, i-1] + a_raw[p,m,i]*tau[m];
      }
      Y_hat[p,m,1:num_years] = a0[m]*to_vector(X) + to_vector(B'*a[p,m, 1:num_knots]);
      P[p, m, 1, 1:num_years] = to_row_vector(inv_logit(Y_hat[p,m,1:num_years])) ;
      P[p, m, 2, 1:num_years] = to_row_vector(1 - P[p, m, 1, 1:num_years]) ;
    }
  }
}

model {
  // Priors
  a0 ~ normal(0, 2);
  tau ~ normal(0, 2);
  sigma_y ~ normal(0, 2);
  
  for(p in 1:P_count){
    for(m in 1:M_count){
      a_raw[p,m] ~ normal(0, 2);
    }
  }

  // Likelihood
  for (k in 1:num_data) {
    Y[k] ~ normal(Y_hat[matchsubnat[k], matchmethod[k], matchyears[k]], sigma_y[matchmethod[k]]);
  }
}

