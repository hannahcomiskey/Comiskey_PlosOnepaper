library("MASS")
set.seed(42)
D <-3
P <- 10 
N <-300

mu_theta <-rep(0,D) # the mean of eta
mu_epsilon<-rep(0,P) # the mean of epsilon
Phi<-diag(rep(1,D))
Psi <- diag(c(0.2079, 0.19, 0.1525, 0.20, 0.36, 0.1875, 0.1875, 1.00, 0.27, 0.27))
l1 <- c(0.99, 0.00, 0.25, 0.00, 0.80, 0.00, 0.50, 0.00, 0.00, 0.00)
l2 <- c(0.00, 0.90, 0.25, 0.40, 0.00, 0.50, 0.00, 0.00, -0.30, -0.30)
l3<-  c(0.00, 0.00, 0.85, 0.80, 0.00, 0.75, 0.75, 0.00, 0.80, 0.80)
L <-cbind(l1,l2,l3) # the loading matrix

Theta <-mvrnorm(N, mu_theta, Phi) # sample factor scores
Epsilon <-mvrnorm(N, mu_epsilon, Psi) # sample error vector
Y<-Theta%*%t(L)+Epsilon# generate observable data

fa.model <- " 
data {
  int<lower=1> N;                // number of 
  int<lower=1> P;                // number of 
  matrix[N,P] Y;                 // data matrix of order [N,P]
  int<lower=1> D;              // number of latent dimensions 
}
transformed data {
  int<lower=1> M;
  vector[P] mu;
  M  <- D*(P-D)+ D*(D-1)/2;  // number of non-zero loadings
  mu <- rep_vector(0.0,P);
}
parameters {    
  vector[M] L_t;   // lower diagonal elements of L
  vector<lower=0>[D] L_d;   // lower diagonal elements of L
  vector<lower=0>[P] psi;         // vector of variances
  real<lower=0>   mu_psi;
  real<lower=0>  sigma_psi;
  real   mu_lt;
  real<lower=0>  sigma_lt;
}
transformed parameters{
  cholesky_factor_cov[P,D] L;  //lower triangular factor loadings Matrix 
  cov_matrix[P] Q;   //Covariance mat
{
  int idx1;
  int idx2;
  real zero; 
  zero <- 0;
  for(i in 1:P){
    for(j in (i+1):D){
      idx1 <- idx1 + 1;
      L[i,j] <- zero; //constrain the upper triangular elements to zero 
    }
  }
  for (j in 1:D) {
      L[j,j] <- L_d[j];
    for (i in (j+1):P) {
      idx2 <- idx2 + 1;
      L[i,j] <- L_t[idx2];
    } 
  }
} 
Q<-L*L'+diag_matrix(psi); 
}
model {
// the hyperpriors 
   mu_psi ~ cauchy(0, 1);
   sigma_psi ~ cauchy(0,1);
   mu_lt ~ cauchy(0, 1);
   sigma_lt ~ cauchy(0,1);
// the priors 
  L_d ~ cauchy(0,3);
  L_t ~ cauchy(mu_lt,sigma_lt);
  psi ~ cauchy(mu_psi,sigma_psi);
//The likelihood
for( j in 1:N)
    Y[j] ~ multi_normal(mu,Q); 
}
"

library("rstan")
library("parallel")

fa.data <-list(P=P,N=N,Y=Y,D=D)

# a function to generate intial values that are slightly jittered for each chain.
init_fun = function() {
  init.values<-list(L_t=rep(0,24)+runif(1,-.1,.1),
                    L_d=rep(.5,D)+runif(1,-.1,.1),
                    psi=rep(.2,P)+runif(1,-.1,.1),
                    sigma_psi=0.15+runif(1,-.1,.1),
                    mu_psi=0.2++runif(1,-.1,.1),
                    sigma_lt=0.5+runif(1,-.1,.1),
                    mu_lt=0.0+runif(1,-.1,.1))
  return(init.values); 
} 

#compile the model
fa.model<- stan(model_code = fa.model, 
                data = fa.data,
                chains =0, 
                pars=c("L","psi","sigma_psi","mu_psi","sigma_lt","mu_lt"))

# run 4 chain in parallel and save samples onto file
# Nchains <- 4
# Niter <- 300
# t_start <- proc.time()[3]
# sflist <- mclapply(1:Nchains, mc.cores = Nchains, 
#                    function(i) stan(fit = fa.model, 
#                                     data =fa.data, 
#                                     pars= c("L","psi","sigma_psi","mu_psi","sigma_lt","mu_lt"), 
#                                     seed = 42,
#                                     iter=Niter,
#                                     init=init_fun,
#                                     #diagnostic_file = paste("diagfile",i,".csv",sep = ""),
#                                     sample_file = paste("sampfile",i,".csv",sep = ""),
#                                     chains = 1, chain_id = i, 
#                                     refresh = -1))
# t_end <- proc.time()[3]

stan(fit = fa.model, 
     data =fa.data, 
     pars= c("L","psi","sigma_psi","mu_psi","sigma_lt","mu_lt"), 
     seed = 42,
     iter=300,
     init=init_fun,
     chains = 1)

t_end <- proc.time()[3]
t_elapsed <- t_end - t_start
fa.fit<- sflist2stanfit(sflist) 
print(fa.fit,probs = c(0.5))
