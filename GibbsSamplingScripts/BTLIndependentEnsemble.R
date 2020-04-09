source(here("GibbsSamplingScripts", "CompositionalBTLGibbs.R"))
set.seed(123)
M1 <- matrix(c(.8, .1, .1,
               .1, .8, .1,
               .1, .1, .1),
             byrow = TRUE, nrow = 3)
M2 <- matrix(c(.5, .25, .25,
               .25, .5, .25,
               .25, .25, .5),
             byrow = TRUE, nrow = 3)
M <- array(NA, dim = c(C, C, 2))
M[,,1] <- M1
M[,,2] <- M2

K <- 2
C <- 3
N_U <- 1000
N_L <- 300
G_L <- sample(1:C, N_L, replace = TRUE, p = rep(1/C, C))
G_U <- sample(1:C, N_U, replace = TRUE, p = c(.25, .5, .25))
A_U <- array(0, dim = c(N_U, C, K))
A_L <- array(0, dim = c(N_L, C, K))
for(k in 1:K) {
    for(i in 1:500) {
        A <- sample(1:C, 1, p = M[G_U[i],,k])
        A_U[i,A,k] <- 1
    }
}
for(k in 1:K) {
    for(i in 1:100) {
        A <- sample(1:C, 1, p = M[G_L[i],,k])
        A_L[i,A,k] <- 1
    }
}
G_L_mat <- matrix(0, nrow = N_L, ncol = C)
for(i in 1:nrow(G_L_mat)){
    G_L_mat[i, G_L[i]] <- 1
}


indiv_samples1 <- compositional_btl(A_U = A_U[,,1],
                                   A_L = A_L[,,1],
                                   G_L = G_L_mat,
                                   causes = as.character(1:C),
                                   thin = 5,
                                   burnin = 1000, ndraws = 10000,
                                   power = 1/100,
                                   epsilon = .01)
indiv_samples2 <- compositional_btl(A_U = A_U[,,2],
                                    A_L = A_L[,,2],
                                    G_L = G_L_mat,
                                    causes = as.character(1:C),
                                    thin = 5,
                                    burnin = 1000, ndraws = 10000,
                                    power = 1/100,
                                    epsilon = .01)
ensemble_samples <- compositional_ensemble_btl(A_U = A_U,
                                               A_L = A_L,
                                               G_L = G_L_mat,
                                               causes = as.character(1:C),
                                               thin = 5,
                                               burnin = 1000, ndraws = 10000,
                                               power = 1/100,
                                               epsilon = .01)
p_1 <- ggs(indiv_samples1$posterior, "p")
p_2 <- ggs(indiv_samples2$posterior, "p")
ensemble_p <- ggs(ensemble_samples$posterior, "p")
p_1 %>%
    group_by(Parameter) %>%
    summarize(m = mean(value))
p_2 %>%
    group_by(Parameter) %>%
    summarize(m = mean(value))
ensemble_p %>%
    group_by(Parameter) %>%
    summarize(m = mean(value))

preds_1 <- max.col(indiv_samples1$predictions, 'first')
preds_2 <- max.col(indiv_samples2$predictions, 'first')
ensemble_preds <- max.col(ensemble_samples$predictions, 'first')
table(G_U, preds_1)
table(G_U, preds_2)
table(G_U, ensemble_preds)
