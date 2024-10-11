data {
  int<lower=1> M_count;       // Number of methods
  int<lower=1> P_count;       // Number of provinces
  int<lower=1> n_years;       // Number of years
  int<lower=1> n_obs;         // Number of observations
  int<lower=1> H;             // Dimension of delta.k
  int<lower=1> K;             // Total number of categories
  int<lower=1> C_count;       // Number of countries
  int<lower=1> R_count;       // Number of regions
  int<lower=1> matchmethod[n_obs];
  int<lower=1> matchsubnat[n_obs];
  int<lower=1> matchyears[n_obs];
  vector[n_obs] se_prop[n_obs]; // Standard error proportions
  vector[M_count] B_ik[P_count, n_years]; // Input data for B_ik
  int<lower=1> kstar[P_count]; // kstar values for provinces
}

parameters {
  // Variance parameters
  real<lower=0> sd_delta[2, M_count];
  vector<lower=-1, upper=1>[10] rho[2]; // Rho parameters
  matrix[2, M_count] alpha_pms; // Alpha parameters
  matrix[2, M_count, P_count] beta_k; // Beta parameters
  vector[H] delta_k[2, M_count, P_count]; // Delta parameters
  vector[2] mu_delta[2]; // Mean parameters
  cholesky_factor_corr[M_count] L_sigma[2]; // Cholesky factor for covariance
}

model {
  // Priors
  for (s in 1:2) {
    for (m in 1:M_count) {
      for (c in 1:C_count) {
        alpha_cms[s, m, c] ~ normal(beta_r[s, m, matchregion[c]], tau_alpha_cms[s]);
      }
      for (r in 1:R_count) {
        beta_r[s, m, r] ~ normal(beta_world[s, m], tau_beta[s]); 
      }
      beta_world[s, m] ~ normal(0, 0.1);
      
      // Variances
      sigma_alpha_pms[s, m] ~ student_t(1, 0, 1);
      tau_alpha_pms[s, m] = 1 / square(sigma_alpha_pms[s, m]);
      sigma_alpha_cms[s, m] ~ student_t(1, 0, 1);
      tau_alpha_cms[s, m] = 1 / square(sigma_alpha_cms[s, m]);
      sigma_beta[s, m] ~ student_t(1, 0, 1);
      tau_beta[s, m] = 1 / square(sigma_beta[s, m]);
    }
  }
  
  // Variance structure
  for (g in 1:2) {
    for (i in 1:M_count) {
      sd_delta[g, i] ~ student_t(1, 0, 1);
      for (j in 1:M_count) {
        if (i != j) {
          sigma_delta[i, j, g] = rho[j + (g - 1) * 4] * sd_delta[g, i] * sd_delta[g, j];
        } else {
          sigma_delta[i, j, g] = square(sd_delta[g, i]);
        }
      }
    }
    
    // Inverse covariance
    inv_sigma_delta[1:M_count, 1:M_count, g] = inverse(sigma_delta[,,g]);
  }

  // Model Estimates
  for (p in 1:P_count) {
    for (m in 1:M_count) {
      for (t in 1:n_years) {
        z[m, p, t] = alpha_pms[1, m, p] + dot_product(B_ik[p, t], beta_k[1, m, p]);
        r[m, p, t] = alpha_pms[2, m, p] + dot_product(B_ik[p, t], beta_k[2, m, p]);
      }
    }
  }

  // Delta sampling
  for (s in 1:2) {
    for (p in 1:P_count) {
      for (j in 1:H) {
        delta_k[s, 1:M_count, p] ~ multi_normal(mu_delta[s], inv_sigma_delta[,,s]);
      }
    }
  }

  // Likelihood
  for (k in 1:n_obs) {
    y[k, 1] ~ normal(z[matchmethod[k], matchsubnat[k], matchyears[k]], se_prop[k, 1]);
    y[k, 2] ~ normal(logit_C[m, p, t], se_prop[k, 2]);
  }
}

generated quantities {
  // Calculate P, Q, U, and logit.CM
  matrix[3, M_count, P_count, n_years] P;
  for (p in 1:P_count) {
    for (m in 1:M_count) {
      for (t in 1:n_years) {
        P[1, m, p, t] = inv_logit(z[m, p, t]);
        Q[m, p, t] = inv_logit(r[m, p, t]);
        U[m, p, t] = 1 - P[1, m, p, t];
        P[2, m, p, t] = Q[m, p, t] * U[m, p, t];
        logit_C[m, p, t] = log(P[2, m, p, t]) - log(U[m, p, t]);
        P[3, m, p, t] = U[m, p, t] - P[2, m, p, t];
      }
    }
  }
}
