library(gtools)
library(here)
library(ggmcmc)
library(coda)
library(tidyverse)
library(stringr)
library(MCMCvis)
library(rstan)
source(here("GibbsSamplingScripts", "CompositionalBTLGibbs.R"))
source(here("GibbsSamplingScripts", "MultiCauseGibbs.R"))
rstan_options(auto_write = TRUE)
options(mc.cores = nchains)
regression_model <- stan_model("DirichletRegression.stan")


set.seed(123)

C <- 5
p_U <- c(.20, .19, .28, .27, .06)
p_U <- p_U / sum(p_U)
p_L <- rep(1/C, C)
N_U <- 2
N_L <- 500

G_L <- rdirichlet(N_L, alpha = 5*p_L)
G_U <- rdirichlet(N_U, alpha = 5*p_U)


M <- matrix(c(.6, .1, .1, .1, .1,
              .1, .6, .1, .1, .1,
              .1, .1, .6, .1, .1,
              .1, .1, .1, .6, .1,
              .1, .1, .1, .1, .6),
            nrow = C, byrow = TRUE)

### Generate A_U and A_L
A_U <- matrix(NA, nrow = N_U, ncol = C)
A_L <- matrix(NA, nrow = N_L, ncol = C)
for(r in 1:N_U) {
    #theta <- runif(1, 5, 10)
    theta <- 5
    A_U[r,] <- rdirichlet(1, alpha = theta * (t(M) %*% G_U[r,])[,1])
}
for(r in 1:N_L) {
    theta <- runif(1, 5, 10)
    #theta <- 5
    A_L[r,] <- rdirichlet(1, alpha = theta * (t(M) %*% G_L[r,])[,1])
}


data <- list(N = N_L, C = C, A = A_L, G = G_L, alpha = 1, beta = 1, epsilon = .01)
dirichlet_samples <- sampling(regression_model, data = data, chains = 1)

comp_multi_btl_samples <- compositional_multi_btl(A_U = A_U,
                                                  A_L = A_L,
                                                  G_L = G_L,
                                                  causes = as.character(1:C),
                                                  thin = 1,
                                                  burnin = 1000, ndraws = 5000,
                                                  alpha = 1, beta = 1,
                                                  power = 1/100,
                                                  epsilon = .01)
dir_M <- ggs(dirichlet_samples, family = "M\\[.*")
MCMCsummary(dirichlet_samples, params = "M")
comp_multi_M <- ggs(comp_multi_btl_samples, family = "M")
ggs_traceplot(comp_multi_M) +
    facet_wrap(~Parameter, nrow = C)
ggs_density(comp_multi_M) +
    facet_wrap(~Parameter, nrow = C, scales = "free_y")
dir_M %>%
    group_by(Parameter) %>%
    summarize(m = mean(value))
comp_multi_M %>%
    group_by(Parameter) %>%
    summarize(m = mean(value))
MCMCsummary(comp_multi_btl_samples, params = "M")
