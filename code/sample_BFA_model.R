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

library("rstan")
library("parallel")

K_choose_2 = (D * (D - 1)) / 2;

fa.data <-list(P=P,N=N,Y=Y,K=D, K_choose_2=K_choose_2)

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

# run 4 chain in parallel and save samples onto file
Niter <- 300
fit <- stan(model_code = 'model/STAN/sample_model.stan', 
            data = fa.data,
            pars= c("L","psi","sigma_psi","mu_psi","sigma_lt","mu_lt"), 
            seed = 42,
            iter=Niter,
            #diagnostic_file = paste("diagfile",i,".csv",sep = "")
            chains = 3)
