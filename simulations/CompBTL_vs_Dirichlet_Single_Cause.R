index <- as.numeric(commandArgs(trailingOnly = TRUE))
output_dir <- "cluster_output_btl_vs_dir_single_cause"
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
library(rstan)
library(MCMCvis)
nchains <- 3
rstan_options(auto_write = TRUE)
options(mc.cores = nchains)
### Should run compileStanModel.R first to avoid recompilation
dirichlet_model <- stan_model("DirichletMixtureSingleCause.stan")
source(here("GibbsSamplingScripts", "CompositionalBTLGibbs.R"))

setting_df <- expand.grid(model = c("dirichlet-mixture",
                                    "overdispered-dirichlet-mixture"),
                          seed_index = 1:500,
                          p_U_index = 1:4)
setting <- setting_df[index,]
set.seed(123)
seeds <- sample(-1e6:1e6, 500, replace = F)
set.seed(seeds[setting$seed_index])
C <- 5
### p_U_list represents all 4 CSMFs from PHMRC
p_U_list <- vector("list", 4)
p_U_list[[1]] <- c(.2, .19, .27, .27, .07)
p_U_list[[2]] <- c(.11, .11, .40, .29, .09)
p_U_list[[3]] <- c(.09, .18, .52, .19, .02)
p_U_list[[4]] <- c(.13, .30, .35, .19, .03)

M <- matrix(c(.65, .35, 0, 0, 0,
              0, .35, .65, 0, 0,
              .1, .1, .6, .1, .1,
              0, 0, 0, .8, .2,
              0, .4, 0, 0, .6),
            nrow = C, byrow = TRUE) 

### M must be all positive
# M[M==0] <- .001
# for(i in 1:C) {
#     M[i,] <- M[i,] / sum(M[i,])
# }

p_U <- p_U_list[[setting$p_U_index]]
### Normalize just in case
p_U <- p_U / sum(p_U)
p_L <- rep(1/C, C)
N_U <- 1000
N_L <- 300


### Generate Gs
tauCluster <- rbinom(N_L + N_U, 1, .7)
smallTau <- runif(N_L + N_U, .1, 1)
largeTau <- runif(N_L + N_U, 10, 20)
tau <- tauCluster*smallTau + (1-tauCluster) * largeTau
G <- matrix(NA, nrow = N_U + N_L, ncol = C)
for(r in 1:nrow(G)) {
    if(r <= N_U) {
        G[r,] <- t(rmultinom(1, 1, prob = p_U))
    } else {
        G[r,] <- t(rmultinom(1, 1, prob = p_L))   
    }
}

### Generate A_U and A_L
A <- matrix(NA, nrow = N_U + N_L, ncol = C)
for(r in 1:nrow(A)) {
    ### True latent variable
    z <- sample(1:C, 1, prob = G[r,])
    if(setting$model == "dirichlet-mixture") {
        M_row <- M[z,]
        which.0 <- M_row == 0
        A[r,which.0] <- 0
        M_pos <- M_row[!which.0]
        A[r,!which.0] <- rdirichlet(1, alpha = 5 * M_pos) 
    } else {
        M_row <- M[z,]
        which.0 <- M_row == 0
        A[r,which.0] <- 0
        M_pos <- M_row[!which.0]
        tauCluster <- rbinom(1, 1, .5)
        smallTau <- runif(1, .1, 1)
        largeTau <- runif(1, 10, 20)
        tau <- tauCluster*smallTau + (1-tauCluster) * largeTau
        A[r,!which.0] <- rdirichlet(1, alpha = tau * M_pos) 
    }
}


removeZeros <- function(X, eps = .0001) {
    X[X==0] <- eps
    return(X/rowSums(X))
}
A_stan <- removeZeros(A)
G_stan <- max.col(G[-(1:N_U),])

data <- list(N_U = N_U, N_L = N_L, C = C,
             A_U = A_stan[1:N_U,],
             A_L = A_stan[-(1:N_U),],
             G_L = G_stan,
             alpha = 1, beta = 1, epsilon = .01, delta = 1)
#iter <- 500
iter <- 6000
warmup <- 1000
start_dir <- Sys.time()
dir_samples <- sampling(dirichlet_model, data = data, chains = nchains, iter = iter,
                        warmup = warmup)
end_dir <- Sys.time()
total_dir <- as.numeric(end_dir-start_dir, units = "secs")
start_comp_btl <- Sys.time()
comp_btl_samples <- compositional_btl(A_U = A[1:N_U,],
                                       A_L = A[-(1:N_U),],
                                       G_L = G[-(1:N_U),],
                                       causes = as.character(1:C),
                                       thin = 1,
                                       burnin = warmup, ndraws = iter,
                                       power = 1/100,
                                       alpha = 1, beta = 1,
                                       epsilon = .01)
end_comp_btl <- Sys.time()
total_comp_btl <- as.numeric(end_comp_btl-start_comp_btl, units = "secs")

true_p <- data.frame(Parameter = paste0("p[", 1:C, "]"),
                     true_value = p_U)
true_M <- data.frame(Parameter = paste0(paste0("M[", rep(1:C, C), ","),
                                        rep(1:C, each = C), "]"),
                     true_value = as.vector(M))
p_dir <- ggs(dir_samples, family = "p\\[.*") %>%
    mutate(method = "Dirichlet")
M_dir <- ggs(dir_samples, family = "M\\[.*") %>%
    mutate(method = "Dirichlet")
p_btl <- ggs(comp_btl_samples, family = "p") %>%
    mutate(method = "Compositional BTL") 
M_btl <-  ggs(comp_btl_samples, family = "M") %>%
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

### Get Rhat for both
dir_rhat <- MCMCsummary(dir_samples, params = c("p", "M")) %>%
    as.data.frame() %>%
    select(Rhat, n.eff) %>%
    rownames_to_column("Parameter") %>%
    mutate(method = "Dirichlet")
btl_rhat <- MCMCsummary(comp_btl_samples, params = c("p","M"))%>%
    as.data.frame() %>%
    select(Rhat, n.eff) %>%
    rownames_to_column("Parameter") %>%
    mutate(method = "Compositional BTL")
rhat_df <- bind_rows(dir_rhat, btl_rhat)
all_results <- inner_join(all_results, rhat_df)
all_results<- cbind(as.data.frame(all_results), setting)

### Finally add time
time_df <- data.frame(method = c("Dirichlet", "Compositional BTL"),
                      time = c(total_dir, total_comp_btl))
all_results <- left_join(all_results, time_df)
saveRDS(all_results, output_file)
quit('no')
