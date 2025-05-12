library("splines")
library("rstan")
source('code/stan_utility.R')

X <- seq(from=-5, to=5, by=.1) # generating inputs
B <- t(bs(X, knots=seq(-5,5,1), degree=3, intercept = TRUE)) # creating the B-splines
num_data <- length(X)
num_basis <- nrow(B)
a0 <- 0.2 # intercept
a <- rnorm(num_basis, 0, 1) # coefficients of B-splines
Y_true <- as.vector(a0*X + a%*%B) # generating the output
Y <- Y_true + rnorm(length(X),0,.2) # adding noise
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
sm<-stan_model("model/STAN/fit_basis.stan")
fit<-sampling(sm,iter=500, data = list(num_data = num_data, num_basis = num_basis, Y = Y, X=X, B=B), control = list(adapt_delta=0.95))

check_all_diagnostics(fit)
