library(gtools)
library(here)
library(ggmcmc)
library(coda)
library(tidyverse)
library(stringr)
source(here("GibbsSamplingScripts", "CompositionalBTLGibbs.R"))
index <- 2
setting_df <- expand.grid(M = paste0("M", 1:3),
                          seed_index = 1:500)
setting <- setting_df[index,]
set.seed(123)
seeds <- sample(-1e6:1e6, 500, replace = F)
set.seed(seeds[setting$seed_index])
C <- 5
p_U <- c(.20, .19, .28, .27, .06)
p_U <- p_U / sum(p_U)
p_L <- rep(1/C, C)
N_U <- 1000
N_L <- 300
G_L <- sample(1:C, N_L, replace = TRUE, p_L)
G_U <- sample(1:C, N_U, replace = TRUE, p_U)


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

thetas <- c(.5, 1, 2, 5, 10)
### Generate A_U and A_L
A_U <- matrix(NA, nrow = N_U, ncol = C)
A_L <- matrix(NA, nrow = N_L, ncol = C)
for(r in 1:N_U) {
    A_U[r,] <- rdirichlet(1, alpha = thetas[G_U[r]] * M[G_U[r],])
}
for(r in 1:N_L) {
    A_L[r,] <- rdirichlet(1, alpha = thetas[G_L[r]] * M[G_L[r],])
}

G_L_mat <- matrix(0, nrow = N_L, ncol = C)
for(i in 1:nrow(G_L_mat)){
    G_L_mat[i, G_L[i]] <- 1
}
### Change matrix for top cause
top_U <- max.col(A_U, 'first')
A_U_top <- matrix(0, nrow = nrow(A_U), ncol = ncol(A_U))
for(i in 1:nrow(A_U_top)) {
    A_U_top[i,top_U[i]] <- 1
}
top_L <- max.col(A_L, 'first')
A_L_top <- matrix(0, nrow = nrow(A_L), ncol = ncol(A_L))
for(i in 1:nrow(A_L_top)) {
    A_L_top[i,top_L[i]] <- 1
}
top_btl_samples <- compositional_btl(A_U = A_U_top,
                                     A_L = A_L_top,
                                     G_L = G_L_mat,
                                     causes = as.character(1:C),
                                     thin = 5,
                                     burnin = 1000, ndraws = 10000,
                                     power = 1/100,
                                     epsilon = .01)
comp_btl_samples <- compositional_btl(A_U = A_U,
                                      A_L = A_L,
                                      G_L = G_L_mat,
                                      causes = as.character(1:C),
                                      thin = 5,
                                      burnin = 1000, ndraws = 10000,
                                      power = 1/100,
                                      epsilon = .01)
uncalib_btl_samples <- compositional_btl(A_U = A_U,
                                         causes = as.character(1:C),
                                         thin = 5,
                                         burnin = 1000, ndraws = 10000,
                                         power = 1/100,
                                         epsilon = .01)
log_lik_gbql_U <- function(A_U, p, M) {
    q_mean <- t(M) %*% p
    return(as.vector(A_U %*% log(q_mean)))
}
log_lik_gbql_L <- function(A_L, G_L, M) {
    q_mean <- G_L %*% M
    return(rowSums(A_L * log(q_mean)))
}


gbql_log_lik <- function(post_samples, A_U, A_L, G_L) {
    log_lik_list <- lapply(1:length(post_samples), function(i) {
        chain_samples <- post_samples[[i]]
        param_names <- colnames(chain_samples)
        p_samples <- chain_samples[,grepl("p", param_names)]
        M_samples <- chain_samples[,grepl("M", param_names)]
        N <- nrow(A_U)
        n <- nrow(A_L)
        log_lik_chain <- matrix(NA, nrow = nrow(chain_samples), ncol = N + n)
        for(s in 1:nrow(log_lik_chain)) {
            p <- p_samples[s,]
            M_vec <- M_samples[s,]
            M_mat <- M_vec
            dim(M_mat) <- c(C, C)
            log_lik_chain[s,1:N] <- log_lik_gbql_U(A_U, p, M_mat)
            log_lik_chain[s,(N+1):(N+n)] <- log_lik_gbql_L(A_L, G_L, M_mat)
        }
        return(log_lik_chain)
    })
    log_lik_mat <- do.call(rbind, log_lik_list)
    return(log_lik_mat)
}

uncalib_log_lik <- function(post_samples, A_U, A_L, G_L, delta = 1, eps = .001) {
    log_lik_list <- lapply(1:length(post_samples), function(i) {
        chain_samples <- post_samples[[i]]
        param_names <- colnames(chain_samples)
        S <- nrow(chain_samples)
        v <- colSums(A_U)
        #p_samples <- chain_samples[,grepl("p", param_names)]
        #M_samples <- chain_samples[,grepl("M", param_names)]
        p_samples <- rdirichlet(S, v + delta)
        C <- ncol(A_U)
        M_mat <- (1-eps) * diag(1, C) + eps / C
        N <- nrow(A_U)
        n <- nrow(A_L)
        log_lik_chain <- matrix(NA, nrow = nrow(chain_samples), ncol = N + n)
        for(s in 1:nrow(log_lik_chain)) {
            p <- p_samples[s,]
            log_lik_chain[s,1:N] <- log_lik_gbql_U(A_U, p, M_mat)
            log_lik_chain[s,(N+1):(N+n)] <- log_lik_gbql_L(A_L, G_L, M_mat)
        }
        return(log_lik_chain)
    })
    log_lik_mat <- do.call(rbind, log_lik_list)
    return(log_lik_mat)
}

log_lik_mat_top <- gbql_log_lik(top_btl_samples,
                                A_U = A_U_top,
                                A_L = A_L_top,
                                G_L = G_L_mat)
log_lik_mat_top_uncalib <- uncalib_log_lik(top_btl_samples,
                                           A_U = A_U_top,
                                           A_L = A_L_top,
                                           G_L = G_L_mat)
log_lik_mat_comp <- gbql_log_lik(comp_btl_samples,
                                A_U = A_U,
                                A_L = A_L,
                                G_L = G_L_mat)
log_lik_mat_comp_uncalib <- uncalib_log_lik(comp_btl_samples,
                                            A_U = A_U,
                                            A_L = A_L,
                                            G_L = G_L_mat)
### All observations
waic(log_lik_mat_top)
waic(log_lik_mat_top_uncalib)

### Just unlabeled set
waic(log_lik_mat_top[,1:nrow(A_U)])
waic(log_lik_mat_top_uncalib[,1:nrow(A_U)])

### Just labeled set
waic(log_lik_mat_top[,-(1:nrow(A_U))])
waic(log_lik_mat_top_uncalib[,-(1:nrow(A_U))])


waic(log_lik_mat_comp[,1:nrow(A_U)])
waic(log_lik_mat_comp_uncalib[,1:nrow(A_U)])

