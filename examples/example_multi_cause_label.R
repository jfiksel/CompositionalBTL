library(here)
library(gtools)
library(here)
library(ggmcmc)
library(coda)
library(tidyverse)
library(stringr)
library(MCMCpack)
source(here("GibbsSamplingScripts", "CompositionalBTLGibbs.R"))
source(here("GibbsSamplingScripts", "MultiCauseGibbs.R"))


id <- 4001
seeds <- 1:500
countries <- c("India", "Mexico", "Philippines", "Tanzania")
setting.df <- expand.grid(seed.index = seeds,
                          country = countries,
                          ncalib = c(25, 100, 200, 400))
setting <- setting.df[id,]
seed.index <- setting$seed.index
c <- as.character(setting$country)


phmrc <- readRDS(here("data", "phmrc_adult.rds"))
country.df <- data.frame(site = c("AP", "Bohol", "Dar", "Mexico", "Pemba", "UP"),
                         country = c("India", 
                                     "Philippines",
                                     "Tanzania",
                                     "Mexico",
                                     "Tanzania",
                                     "India"))
country <- country.df$country[match(phmrc$site, country.df$site)]
test <- phmrc[country == c,]
external <- c("Road Traffic",
              "Falls",
              "Homicide",
              "Suicide",
              "Fires",
              "Drowning",
              "Other Injuries",
              "Poisonings",
              "Bite of Venomous Animal")
circulatory <- c("Stroke",
                 "Other Cardiovascular Diseases",
                 "Acute Myocardial Infarction")
non_communicable <- c("Other Non-communicable Diseases",
                      "Colorectal Cancer",
                      "Breast Cancer",
                      "Leukemia/Lymphomas",
                      "Prostate Cancer",
                      "Esophageal Cancer",
                      "Stomach Cancer",
                      "Lung Cancer",
                      "Cervical Cancer",
                      "Renal Failure",
                      "Epilepsy",
                      "Cirrhosis",
                      "COPD",
                      "Diabetes",
                      "Asthma")
infectious <- c("Malaria",
                "Pneumonia",
                "Diarrhea/Dysentery",
                "AIDS",
                "TB",
                "Other Infectious Diseases")
maternal <- c("Maternal")

cause_map <- data.frame(causes = c(external, circulatory, non_communicable,
                                   infectious, maternal),
                        broad_cause = rep(c("external", "circulatory",
                                            "non_communicable", "infectious", "maternal"),
                                          c(length(external), length(circulatory),
                                            length(non_communicable), length(infectious),
                                            length(maternal))))
causes <- c("external", "circulatory","non_communicable", "infectious", "maternal")
C <- length(causes)


### Get model probabilities
models <- c("insilico", "tariff", "interva", "nbc")
models_list <- lapply(models, function(model) {
    model_file <- here("train_models", "adult_models", paste0(model, "_model_", c, "_probs.rds"))
    model_probs <- readRDS(model_file)
    
    model_broad_probs <- model_probs
    colnames(model_broad_probs) <- cause_map$broad_cause[match(colnames(model_broad_probs),
                                                               cause_map$causes)]
    
    mn <- model.matrix(~ colnames(model_broad_probs) + 0)
    model_broad_probs <- model_broad_probs %*% mn
    colnames(model_broad_probs) <- gsub("colnames\\(model_broad_probs\\)", "",
                                        colnames(model_broad_probs))
    model_broad_probs <- model_broad_probs[,causes]
    
    model_top_cause <- colnames(model_probs)[max.col(model_probs, 'first')]
    model_top_cause <- cause_map$broad_cause[match(model_top_cause, cause_map$causes)]
    top_cause_mat <- matrix(0, nrow = nrow(model_probs), ncol = C)
    for(i in 1:nrow(top_cause_mat)) {
        cause_index <- which(causes == model_top_cause[i])
        top_cause_mat[i,cause_index] <- 1
    }
    
    model_top_broad_cause <- max.col(model_broad_probs, 'first')
    top_broad_cause_mat <- matrix(0, nrow = nrow(model_broad_probs), ncol = C)
    for(i in 1:nrow(top_broad_cause_mat)) {
        top_broad_cause_mat[i,model_top_broad_cause[i]] <- 1
    }
    return(list(comp_mat = model_broad_probs,
                top_mat = top_cause_mat,
                top_broad_mat = top_broad_cause_mat))
})

### Create array
K <- length(models_list)
comp_array <- array(dim = c(nrow(test), C, K))
top_array <- array(dim = c(nrow(test), C, K))
top_broad_array <- array(dim = c(nrow(test), C, K))
for(k in 1:K) {
    comp_array[,,k] <- models_list[[k]]$comp_mat
    top_array[,,k] <- models_list[[k]]$top_mat
    top_broad_array[,,k] <- models_list[[k]]$top_broad_mat
}

true_broad_cause <- cause_map$broad_cause[match(test$gs_text34, cause_map$causes)]
true_p <- sapply(causes, function(c) mean(true_broad_cause == c))

true_df <- data.frame(true_value = true_p,
                      Parameter = paste0("p[", 1:C, "]"),
                      cause_name = causes)

#########################################
### Create pairs and uncertainty G matrix
set.seed(123)
all_obs <- 1:length(true_broad_cause)
obs1 <- c()
obs2 <- c()
is_paired <- rep(FALSE, length(all_obs))
i <- 1
resample <- function(x, ...) x[sample.int(length(x), ...)]
while(sum(!is_paired) > 0) {
    pair <- resample(all_obs[!is_paired], 2, replace = TRUE)
    obs1 <- c(obs1, pair[1])
    obs2 <- c(obs2, pair[2])
    is_paired[pair] <- TRUE
    i <- i + 1
}

true_broad_cause <- as.character(true_broad_cause)
G <- matrix(0, nrow = length(all_obs), ncol = C)
for(i in 1:length(obs1)) {
    cause1 <- which(causes == true_broad_cause[obs1[i]])
    cause2 <- which(causes == true_broad_cause[obs2[i]])
    G[c(obs1[i], obs2[i]), cause1] <- G[c(obs1[i], obs2[i]), cause1] + .5
    G[c(obs1[i], obs2[i]), cause2] <- G[c(obs1[i], obs2[i]), cause2] + .5
}


set.seed(123)
seeds <- sample(-1e6:1e6, size = 500, replace = F)
init.seed = seeds[seed.index]
set.seed(init.seed)
### Sample using inverse probabilities
sample_probs <- 1 / true_df[match(true_broad_cause, true_df$cause_name), "true_value"]
calib_indices <- sample(1:nrow(test), setting$ncalib, prob = sample_probs, replace = FALSE)
G_L <- G[calib_indices,]

csmf_acc <- function(truth, csmf) {
    acc <- 1 - sum(abs(truth - csmf))/2/(1 - min(truth))
    return(acc)
}

ndraw_gbql <- 10000

### Get draws for insilico
k <- 1
insilico_comp_btl <- compositional_multi_btl(A_U = comp_array[,,k],
                                             A_L = comp_array[calib_indices,,k],
                                             G_L = G_L,
                                             causes = as.character(1:C),
                                             thin = 1,
                                             burnin = 1000, ndraws = ndraw_gbql,
                                             power = 1/100,
                                             epsilon = .01,
                                             alpha = 10)
### Get calibrated predictions
insilico_calib_p <-
    ggs(insilico_comp_btl, "p") %>%
    group_by(Parameter) %>%
    summarise(estimate = mean(value)) 

### Ensemble

ens_comp_btl <- compositional_multi_ensemble_btl(A_U = comp_array,
                                                 A_L = comp_array[calib_indices,,],
                                                 G_L = G_L,
                                                 causes = as.character(1:C),
                                                 thin = 1,
                                                 burnin = 1000, ndraws = ndraw_gbql,
                                                 power = 1/100,
                                                 epsilon = .01,
                                                 alpha = 10)

ens_top_comp <-
    ggs(ens_comp_btl, "p") %>%
    group_by(Parameter) %>%
    summarise(estimate = mean(value)) 
