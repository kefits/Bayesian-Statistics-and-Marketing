library(rstan)
library(dplyr)
library(bayesm)
library(mvtnorm)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

addSelectFlag <- function(x) {
  len_x <- length(x)
  x_max <- max(x)
  
  res <- rep(0, len_x+1)
  if (sum(x < 0) == len_x) {
    res[len_x+1] <- 1
  } else {
    for(i in 1:len_x) {
      if (x[i] == x_max) {
        res[i] <- 1
      } else {
        res[i] <- 0
      }
    }
  }
  return(res)
}

N = 1600
p = 6
beta = rep(2, p-1)
sigma = sqrt(c(1, 2, 3, 4, 5))
iota = rep(1, p-1)
rho = 0.5
Sigma = diag(sigma) %*% (rho * iota %*% t(iota) + (1-rho) * diag(rep(1, p-1))) %*% diag(sigma)

X <- matrix(runif(N*(p-1)^2, -2, 2), N*(p-1), p-1)
for(i in 1:N) {
  X[i,] <- X[i,] - runif(1, -2, 2)
}

mu <- X %*% beta

y <- matrix(rep(0, p*N), N, p)
for(i in 1:N){
  w <- rmvnorm(1, mu[((p-1)*(i-1)+1):((p-1)*i)], Sigma)
  y[i,] <- addSelectFlag(w)
}


params <- c("beta_tilde", "w", "Sigma_tilde")
A=diag(rep(100, dim(X)[2]))
nu=8
d.data <- list(N=N, M=dim(X)[2], D=(p-1), Y=y, X=X,
               nu=nu, A=A)
d.model <- stan_model("MNP.stan")
d.fit <- sampling(d.model, data=d.data, iter = 100, chains = 4, pars = params)

print(d.fit, "beta_tilde")

draws <- extract(d.fit)
w <- draws$w %>% 
  as.data.frame()
beta_tilde <- draws$beta_tilde %>% 
  as.data.frame()
Sigma_tilde <- draws$Sigma_tilde %>% 
  as.data.frame()

cor(w$`1.4`, w$`1.5`)
plot(beta_tilde$V1, type="l")
plot(Sigma_tilde$`2.1`, type="l")


#############################################################


N = 1600
p = 6
beta = rep(2,1) 
sigma = sqrt(c(1, 2, 3, 4, 5))
iota = rep(1, p-1)
rho = 0.5
Sigma = diag(sigma) %*% (rho * iota %*% t(iota) + (1-rho) * diag(rep(1, p-1))) %*% diag(sigma)


X <- matrix(runif(N*p, -2, 2), N, p)
X <- createX(p=p, na=1, nd=NULL, Xa=X, Xd=NULL, INT = F, DIFF = T)

mu <- X %*% beta

y <- matrix(rep(0, p*N), N, p)
for(i in 1:N){
  w <- rmvnorm(1, mu[((p-1)*(i-1)+1):((p-1)*i)], Sigma)
  y[i,] <- addSelectFlag(w)
}
y <- y %*% seq(1,6)
R <- 2000
Data <- list(p=p,X=X,y=y)
Mcmc <- list(R=R)
out <- rmnpGibbs(Data=Data, Mcmc=Mcmc)
Sigma_bm <- as.data.frame(out$sigmadraw[1:R,])



