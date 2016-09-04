library(bayesm)
library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data("margarine")

hhid_selected <- margarine$choicePrice %>% 
  filter(choice %in% c(1,2,3,4,5,7)) %>% 
  group_by(hhid) %>% 
  summarise(purc_cnt = n()) %>% 
  filter(purc_cnt >= 5)

choicePrice.selected <- margarine$choicePrice %>% 
  filter(choice %in% c(1,2,3,4,5,7) & hhid %in% hhid_selected$hhid)
choicePrice.selected$choice[choicePrice.selected$choice == 7] <- 6

y <- choicePrice.selected$choice

X <- choicePrice.selected %>% 
  select(3,4,5,6,7,9)


d.dat <- list(N_x=nrow(X), p_x=ncol(X), p_y=6,
              y=y, X=X, K=5)
d.fit <- stan("MNL_mixture.stan", data = d.dat, iter = 500, chains = 4)

print(d.fit, pars = "pi",
      probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
      digits_summary = 2, include = TRUE)



draws <- extract(d.fit)
beta <- as.data.frame(draws$beta)
V_b <- as.data.frame(draws$V_b)


# Plot
plotDensities <- function(varname){
  par(mfrow = c(2,3))
  for(i in 1:6){
    plot(density(varname[[i]]))
  }
  par(mfrow = c(1,1))
}

plotDensities(beta)

