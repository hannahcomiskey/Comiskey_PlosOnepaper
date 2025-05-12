data {
  int num_data;
  int num_basis;
  int num_years;
  vector[num_data] Y;
  vector[num_years] X;
  matrix[num_basis, num_years] B;
  int matchmethod[num_data] ; // method indexing
  int matchyears[num_data]; // year indexing
  int matchsubnat[num_data]; // subnat indexing 
  int M_count; // Number of methods
  int P_count; // Number of provinces
  int S_count; // Number of sectors
}
parameters {
  row_vector[num_basis] a_raw[P_count, M_count];
  vector[M_count] a0[P_count];
  vector<lower=0>[M_count] sigma_y;
  vector<lower=0>[M_count] tau;
}
transformed parameters {
  row_vector[num_basis] a[P_count, M_count];
  vector[num_years] Y_hat[P_count, M_count];
  matrix<lower=0, upper=1>[S_count, num_years] P[P_count, M_count]; // logit observation
  
  for(p in 1:P_count) {
    for(m in 1:M_count) { 
      a[p,m,1] = a_raw[p,m,1]*tau[m];
      for (i in 2:(num_basis-1)) {
        a[p,m,i] = a[p,m, i-1] + a_raw[p,m,i]*tau[m]; // penalising splines
      }
      a[p,m,num_basis] = -sum(a[p,m,1:(num_basis-1)]); # hard sum to 0 constraint
      Y_hat[p,m,1:num_years] = a0[p,m]*X + to_vector(a[p,m, 1:num_basis]*B);
      P[p, m, 1, 1:num_years] = to_row_vector(inv_logit(Y_hat[p,m,1:num_years])) ;
      P[p, m, 2, 1:num_years] = to_row_vector(1 - P[p, m, 1, 1:num_years]) ;
    }
  }
}
model {
  for(p in 1:P_count) {
    a0[p] ~ cauchy(0, 1);
    for(m in 1:M_count){ 
      a_raw[p, m] ~ normal(0, 1);
    }
  }
  tau ~ cauchy(0, 1);
  sigma_y ~ cauchy(0, 1);
  
  
  // Likelihood
  for (k in 1:num_data) {
    Y[k] ~ normal(Y_hat[matchsubnat[k], matchmethod[k], matchyears[k]], sigma_y[matchmethod[k]]);
  }
}


