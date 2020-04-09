index <- as.numeric(commandArgs(trailingOnly = TRUE))
output_dir <- "cluster_output_transerror"
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
library(openVA)
source(here("GibbsSamplingScripts", "CompositionalBTLGibbs.R"))
setting_df <- expand.grid(M = paste0("M", 1:3),
                          seed_index = 1:500)
setting <- setting_df[index,]
set.seed(123)
seeds <- sample(-1e6:1e6, 500, replace = F)

### Generate conditional probability symptoms and source data
### rows are symptoms, columns are causes
set.seed(123)
C <- 5
N_S <- 5000
G_Source <- sample(1:C, N_S, replace = TRUE, c(.11, .18, .43, .23, .05))
prob_mat_source <- matrix(NA, nrow = 21, ncol = C)
for(i in 1:C) {
    prob_mat_source[,i] <- runif(21)
}
train_mat <- matrix(NA, nrow = N_S, ncol = 21)
for(i in 1:nrow(train_mat)){
    train_mat[i,] <- rbinom(21, size = 1, prob = prob_mat_source[,G_Source[i]])
}
train <- data.frame(train_mat)
train[train==1] <- "Y"
train[train==0] <- ""
train <- cbind(deathID = as.character(1:N_S), cause = as.character(G_Source), train)


if(setting$M == "M1") {
    prob_mat_target <- prob_mat_source
} else if (setting$M == "M2") {
    ### Induce transfer error similar to M2 in previous simulation
    prob_mat_target <- prob_mat_source
    prob_mat_target[5:14,2] <- prob_mat_source[5:14,1]
    prob_mat_target[1:2,3] <- prob_mat_source[1:2, 1]
    prob_mat_target[3:4,3] <- prob_mat_source[3:4, 2]
    prob_mat_target[5:6,3] <- prob_mat_source[5:6, 4]
    prob_mat_target[7:8,3] <- prob_mat_source[7:8, 5]
    prob_mat_target[12:20,4] <- prob_mat_source[12:20,3]
}else {
    ### Induce transfer error similar to M3 (small misclassifications)
    prob_mat_target <- prob_mat_source
    for(i in 1:C) {
        diffvec <- setdiff(1:C, i)
        for(j in seq_along(diffvec)) {
            prob_mat_target[(j*2-1):(j*2), i] <- prob_mat_source[(j*2-1):(j*2), diffvec[j]] 
        }
    }
}


set.seed(seeds[setting$seed_index])
p_U <- c(.20, .19, .28, .27, .06)
p_U <- p_U / sum(p_U)
p_L <- rep(1/C, C)
N_U <- 1000
N_L <- 300
G_L <- sample(1:C, N_L, replace = TRUE, p_L)
G_U <- sample(1:C, N_U, replace = TRUE, p_U)


testmat_L <- matrix(NA, nrow = N_L, ncol = 21)
testmat_U <- matrix(NA, nrow = N_U, ncol = 21)

for(i in 1:nrow(testmat_U)){
    testmat_U[i,] <- rbinom(21, size = 1, prob = prob_mat_target[,G_U[i]])
}
for(i in 1:nrow(testmat_L)){
    testmat_L[i,] <- rbinom(21, size = 1, prob = prob_mat_target[,G_L[i]])
}
testmat <- rbind(testmat_U, testmat_L)
testmat[testmat==1] <- "Y"
testmat[testmat==0] <- ""
test <- data.frame(testmat)
test <- cbind(deathID = as.character(1:(N_L + N_U)),
              test)

insilicova_U <- codeVA(data = test[1:N_U,], data.type = "customize",
                     model = "InSilicoVA",
                     data.train = train, causes.train = "cause",
                     jump.scale = 0.05,
                     Nsim=5000, auto.length = FALSE)
insilicova_L <- codeVA(data = test[-(1:N_U),], data.type = "customize",
                     model = "InSilicoVA",
                     data.train = train, causes.train = "cause",
                     jump.scale = 0.05,
                     Nsim=5000, auto.length = FALSE)
A_U <- getIndivProb(insilicova_U)
A_U <- A_U[,as.character(1:5)]
A_L <- getIndivProb(insilicova_L)
A_L <- A_L[,as.character(1:5)]

### Create G_L_mat
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

### Compare estimates of p
csmf_acc <- function(truth, csmf) {
    acc <- 1 - sum(abs(truth - csmf))/2/(1 - min(truth))
    return(acc)
}
true_p <- data.frame(Parameter = paste0("p[", 1:5, "]"),
                     truth = p_U)
raw_p <- data.frame(Parameter = paste0("p[", 1:5, "]"),
                    raw_mean = colMeans(A_U))
top_btl_p <- ggs(top_btl_samples$posterior, family = "^p\\[.*") %>%
    group_by(Parameter) %>%
    summarize(top_post_mean = mean(value))
comp_btl_p <- ggs(comp_btl_samples$posterior, family = "^p\\[.*") %>%
    group_by(Parameter) %>%
    summarize(comp_post_mean = mean(value))

p_df <- inner_join(comp_btl_p, top_btl_p) %>%
    inner_join(true_p) %>%
    inner_join(raw_p)

top_btl_csmf_acc = csmf_acc(p_df$truth, p_df$top_post_mean)
raw_csmf_acc = csmf_acc(p_df$truth, p_df$raw_mean)
comp_btl_csmf_acc = csmf_acc(p_df$truth, p_df$comp_post_mean)

### individual predictions
getCCC <- function(fitted, truth){
    fitted <- as.character(fitted) 
    truth <- as.character(truth)
    C <- length(unique(truth))
    cccj <- rep(NA, C)
    correct <- fitted[which(fitted == truth)]
    N <- length(truth)
    for(i in 1:length(unique(truth))){
        c <- sort(unique(truth))[i]
        if(length(which(truth == c)) == 0) next
        cccj[i] <- length(which(correct == c))/length(which(truth == c))
        cccj[i] <- (cccj[i] - 1/C) / (1 - 1/C)
    }
    ccc <- mean(cccj, na.rm = TRUE)
    return(ccc)
}
### get top classes for compositional BTL
top_btl_preds <- top_btl_samples$predictions
comp_btl_preds <- comp_btl_samples$predictions
top_btl_top <- max.col(top_btl_preds, 'first')
top_btl_ccc <- getCCC(top_btl_top, G_U)
comp_btl_top <- max.col(comp_btl_preds, 'first')
comp_btl_ccc <- getCCC(comp_btl_top, G_U)
raw_top <- max.col(A_U, 'first')
raw_ccc <- getCCC(raw_top, G_U)

results_df <- data.frame(acc = c(top_btl_csmf_acc, comp_btl_csmf_acc, raw_csmf_acc,
                                 top_btl_ccc, comp_btl_ccc, raw_ccc),
                         method = rep(c('BTL_Top', 'BTL_Comp', 'Raw'), 2),
                         metric = rep(c('CSMF', 'CCC'), each = 3))
results_df <- cbind(results_df, setting)
saveRDS(results_df, output_file)
quit('no')