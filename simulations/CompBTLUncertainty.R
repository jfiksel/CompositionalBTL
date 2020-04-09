index <- as.numeric(commandArgs(trailingOnly = TRUE))
output_dir <- "cluster_output_uncertanity_vs_regression"
if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
}
output_file <- file.path(output_dir, paste0("run-", index, ".rds"))
if(file.exists(output_file)){
    quit('no')
}
library(gtools)
library(here)
library(ggmcmc)
library(coda)
library(tidyverse)
library(stringr)
source(here("GibbsSamplingScripts", "CompositionalBTLGibbs.R"))
source(here("GibbsSamplingScripts", "MultiCauseGibbs.R"))
setting_df <- expand.grid(M = paste0("M", 1:3),
                          N_L = c(25, 100, 300),
                          seed_index = 1:500)
setting <- setting_df[index,]
set.seed(123)
seeds <- sample(-1e6:1e6, 500, replace = F)
set.seed(seeds[setting$seed_index])
C <- 5
### p_U same as Tanzania
p_U <- c(.11, .11, .40, .29, .09)
p_U <- p_U / sum(p_U)
p_L <- rep(1/C, C)
N_U <- 1000
N_L <- setting$N_L

### Sample G_L and G_U
G_L <- matrix(NA, nrow = N_L, ncol = C)
G_U <- matrix(NA, nrow = N_U, ncol = C)
### Probability of getting a completely certain response for G
pi_L <- rbeta(N_L, 1, 1)
pi_U <- rbeta(N_U, 1, 1)
for(r in 1:nrow(G_L)) {
    isMult <- rbinom(1, 1, prob = pi_L[r])
    if(isMult) {
        G_L[r,] <- rmultinom(1, 1, prob = p_L) 
    } else {
        G_L[r,] <- rdirichlet(1, alpha = 5 * p_L)
    }
}
for(r in 1:nrow(G_U)) {
    isMult <- rbinom(1, 1, prob = pi_U[r])
    if(isMult) {
        G_U[r,] <- rmultinom(1, 1, prob = p_U) 
    } else {
        G_U[r,] <- rdirichlet(1, alpha = 5 * p_U)
    }
}

### Generate true latent states
z_L <- sapply(1:nrow(G_L), function(r) {
    sample(1:C, 1, prob = G_L[r,])
})
z_U <- sapply(1:nrow(G_U), function(r) {
    sample(1:C, 1, prob = G_U[r,])
})


if(setting$M == "M1") {
    M <- diag(1, C)
} else if (setting$M == "M2") {
    M <- matrix(c(1, 0, 0, 0, 0,
                  .65, .35, 0, 0, 0,
                  .1, .1, .6, .1, .1,
                  0, 0, .5, .5, 0,
                  0, 0, 0, 0, 1),
                nrow = C, byrow = TRUE) 
} else {
    M <- matrix(c(.6, .1, .1, .1, .1,
                  .1, .6, .1, .1, .1,
                  .1, .1, .6, .1, .1,
                  .1, .1, .1, .6, .1,
                  .1, .1, .1, .1, .6),
                nrow = C, byrow = TRUE)
} 
### Make it so M cannot be 0 (presents problems for dirichlet)
M[M == 0] <- .001
for(i in 1:nrow(M)) {
    M[i,] <- M[i,] / sum(M[i,])
}

### Generate A_U and A_L
A_U_unc <- A_U_reg <- matrix(NA, nrow = N_U, ncol = C)
A_L_unc <- A_L_reg <- matrix(NA, nrow = N_L, ncol = C)

for(r in 1:N_U) {
    theta <- runif(1, 1, 10)
    A_U_unc[r,] <- rdirichlet(1, alpha = theta * M[z_U[r],])
    A_U_reg[r,] <-  rdirichlet(1, alpha = theta * (t(M) %*% G_U[r,])[,1])
}
for(r in 1:N_L) {
    theta <- runif(1, 1, 10)
    A_L_unc[r,] <- rdirichlet(1, alpha = theta * M[z_L[r],])
    A_L_reg[r,] <-  rdirichlet(1, alpha = theta * (t(M) %*% G_L[r,])[,1])
}


unc_btl_samples <- compositional_multi_btl(A_U = A_U_unc,
                                           A_L = A_L_unc,
                                           G_L = G_L,
                                           causes = as.character(1:C),
                                           thin = 5,
                                           alpha = 15, beta = .5,
                                           burnin = 1000, ndraws = 5000,
                                           power = 1/100,
                                           epsilon = .01)
reg_btl_samples <- compositional_multi_btl(A_U = A_U_reg,
                                           A_L = A_L_reg,
                                           G_L = G_L,
                                           causes = as.character(1:C),
                                           thin = 5,
                                           alpha = 15, beta = .5,
                                           burnin = 1000, ndraws = 5000,
                                           power = 1/100,
                                           epsilon = .01)


true_p <- data.frame(Parameter = paste0("p[", 1:C, "]"),
                     true_value = p_U)
true_M <- data.frame(Parameter = paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), "]"),
                     true_value = as.vector(M))
p_dir <- ggs(dir_samples, family = "p\\[.*") %>%
    mutate(method = "Dirichlet")
M_dir <- ggs(dir_samples, family = "M\\[.*") %>%
    mutate(method = "Dirichlet")
p_btl <- ggs(comp_btl_samples$posterior, family = "p") %>%
    mutate(method = "Compositional BTL") 
M_btl <-  ggs(comp_btl_samples$posterior, family = "M") %>%
    mutate(method = "Compositional BTL")  

p_all <-
    bind_rows(p_dir, p_btl) %>%
    group_by(Parameter, method) %>%
    summarise(mean = mean(value), ci_L = quantile(value, .025), ci_U = quantile(value, .975)) %>%
    inner_join(true_p)
M_all <-
    bind_rows(M_dir, M_btl) %>%
    group_by(Parameter, method) %>%
    summarise(mean = mean(value), ci_L = quantile(value, .025), ci_U = quantile(value, .975)) %>%
    inner_join(true_M)
all_results <- bind_rows(p_all, M_all)

true_p <- data.frame(Parameter = paste0("p[", 1:5, "]"),
                     truth = p_U)
raw_p <- data.frame(Parameter = paste0("p[", 1:5, "]"),
                    raw_mean = colMeans(A_U))
unc_btl_p <-
    ggs(unc_btl_samples, family = "^p\\[.*") %>%
    group_by(Parameter) %>%
    summarize(unc_post_mean = mean(value))

reg_btl_p <-
    ggs(reg_btl_samples, family = "^p\\[.*") %>%
    group_by(Parameter) %>%
    summarize(reg_post_mean = mean(value))

p_df <-
    inner_join(unc_btl_p, reg_btl_p) %>%
    inner_join(true_p) %>%
    inner_join(raw_p)

p_df <- cbind(p_df, setting)

true_M <- data.frame(Parameter = paste0(paste0("M[", rep(1:C, C), ","), rep(1:C, each = C), "]"),
                     true_value = as.vector(M))

unc_btl_p <-
    ggs(unc_btl_samples, family = "^p\\[.*") %>%
    group_by(Parameter) %>%
    summarize(unc_post_mean = mean(value))

reg_btl_p <-
    ggs(reg_btl_samples, family = "^p\\[.*") %>%
    group_by(Parameter) %>%
    summarize(reg_post_mean = mean(value))

saveRDS(p_df, output_file)
quit('no')

