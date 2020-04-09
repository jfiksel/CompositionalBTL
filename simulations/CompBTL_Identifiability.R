library(gtools)
library(here)
library(ggmcmc)
library(coda)
library(tidyverse)
library(stringr)
source(here("GibbsSamplingScripts", "CompositionalBTLGibbs.R"))
set.seed(123)
C <- 5
p_U <- c(.1, .4, .2, .2, .1)
p_U <- p_U / sum(p_U)
p_L <- rep(1/C, C)
N_U <- 2000
N_L <- 300
G_L <- sample(1:C, N_L, replace = TRUE, p_L)
G_U <- sample(1:C, N_U, replace = TRUE, p_U)


### Generate A_U and A_L
A_U <- matrix(0, nrow = N_U, ncol = C)
A_L <- matrix(0, nrow = N_L, ncol = C)
for(r in 1:N_U) {
    if(G_U[r] == 1) {
        A_U[r,1] <- 1
    } else {
        A_U[r,1] <- .7
        A_U[r,G_U[r]] <- .3
    }
}
for(r in 1:N_L) {
    if(G_L[r] == 1) {
        A_L[r,1] <- 1
    } else {
        A_L[r,1] <- .7
        A_L[r,G_L[r]] <- .3
    }
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

top_p <- ggs(top_btl_samples$posterior, family = "p")
comp_p <- ggs(comp_btl_samples$posterior, family = "p")
true_p <- data.frame(Parameter = paste0("p[", 1:5, "]"),
                     truth = p_U)
p_df <-  bind_rows(mutate(top_p, method = "Top Cause"),
                   mutate(comp_p, method = "Compositional")) %>%
    inner_join(true_p)
dens_fig <-
    ggplot(p_df, aes(x = value, color = method))+
    geom_line(stat = 'density') +
    geom_vline(aes(xintercept = truth)) +
    facet_wrap(~Parameter, scales = "free_y")
ggsave(here("simulations", "figs", "density_identifiability.pdf"),
       dens_fig, width = 8, height = 8)

### M samples
top_m <- ggs(top_btl_samples$posterior, family = "M")
ggs_traceplot(top_m) +
    facet_wrap(~Parameter)

comp_m <- ggs(comp_btl_samples$posterior, family = "M")
ggs_traceplot(comp_m) +
    facet_wrap(~Parameter)

