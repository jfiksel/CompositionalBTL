library(here)
cluster_dir <- here("phmrc_adult_5_causes_analysis", "cluster_output")
if(!dir.exists(cluster_dir)){
    dir.create(cluster_dir)
}
id <- as.numeric(commandArgs(trailingOnly = TRUE))
output_file <- file.path(cluster_dir, paste0("run-", id, ".rds"))
if(file.exists(output_file)){
    quit('no')
}

library(gtools)
library(here)
library(ggmcmc)
library(coda)
library(tidyverse)
library(stringr)
library(MCMCpack)
source(here("GibbsSamplingScripts", "CompositionalBTLGibbs.R"))


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
set.seed(123)
seeds <- sample(-1e6:1e6, size = 500, replace = F)
init.seed = seeds[seed.index]
set.seed(init.seed)
### Sample using inverse probabilities
sample_probs <- 1 / true_df[match(true_broad_cause, true_df$cause_name), "true_value"]
calib_indices <- sample(1:nrow(test), setting$ncalib, prob = sample_probs, replace = FALSE)
G_calib = as.character(true_broad_cause[calib_indices])
G_indices <- sapply(G_calib, function(c) which(causes == c))
G_L <- matrix(0, nrow = length(G_indices), ncol = C)
for(i in 1:nrow(G_L)) {
    G_L[i, G_indices[i]] <- 1
}

csmf_acc <- function(truth, csmf) {
    acc <- 1 - sum(abs(truth - csmf))/2/(1 - min(truth))
    return(acc)
}

ndraw_gbql <- 10000
posteriors_list <- lapply(1:length(models), function(k) {
    model <- models[k]
    top_btl <- compositional_btl(A_U = top_array[,,k],
                                 A_L = top_array[calib_indices,,k],
                                 G_L = G_L,
                                 causes = as.character(1:C),
                                 thin = 5,
                                 burnin = 1000, ndraws = ndraw_gbql,
                                 power = 1/100,
                                 epsilon = .01,
                                 alpha = 10)
    comp_btl <- compositional_btl(A_U = comp_array[,,k],
                                       A_L = comp_array[calib_indices,,k],
                                       G_L = G_L,
                                       causes = as.character(1:C),
                                       thin = 5,
                                       burnin = 1000, ndraws = ndraw_gbql,
                                       power = 1/100,
                                       epsilon = .01,
                                  alpha = 10)
    
    ### Get raw predictions
    top_raw <- sapply(1:C, function(c) mean(max.col(top_array[,,k], 'first') == c))
    comp_raw <- sapply(1:C, function(c) mean((comp_array[,c,k])))
    
    raw_df <- data.frame(Parameter = rep(paste0("p[", 1:5, "]"), 2),
                         estimate = c(top_raw, comp_raw),
                         method = rep(c("Top", "Compositional"), each = C),
                         calibrated = FALSE)
    
    ### Get calibrated predictions
    top_calib <-
        ggs(top_btl, "p") %>%
        group_by(Parameter) %>%
        summarise(estimate = mean(value)) %>%
        mutate(method = "Top", calibrated = TRUE)
    top_comp <-
        ggs(comp_btl, "p") %>%
        group_by(Parameter) %>%
        summarise(estimate = mean(value)) %>%
        mutate(method = "Compositional", calibrated = TRUE)
    
    csmf_estimates <-
        raw_df %>%
        bind_rows(top_calib) %>%
        bind_rows(top_comp)
    
    csmf_estimates <-
        inner_join(csmf_estimates, true_df, by = "Parameter") %>%
        mutate(model = models[k])
    csmf_acc_df <-
        csmf_estimates %>%
        group_by(method, calibrated, model) %>%
        summarise(csmfa = csmf_acc(true_value, estimate)) %>%
        mutate(cccsmfa = (csmfa - .632)/(1-.632))
    return(csmf_acc_df)
})

### Now ensemble
ens_top_btl <- compositional_ensemble_btl(A_U = top_array,
                                          A_L = top_array[calib_indices,,],
                                          G_L = G_L,
                                          causes = as.character(1:C),
                                          thin = 5,
                                          burnin = 1000, ndraws = ndraw_gbql,
                                          power = 1/100,
                                          epsilon = .01,
                                          alpha = 10)

ens_comp_btl <- compositional_ensemble_btl(A_U = comp_array,
                                           A_L = comp_array[calib_indices,,],
                                           G_L = G_L,
                                           causes = as.character(1:C),
                                           thin = 5,
                                           burnin = 1000, ndraws = ndraw_gbql,
                                           power = 1/100,
                                           epsilon = .01,
                                           alpha = 10)
ens_top_calib <-
    ggs(ens_top_btl, "p") %>%
    group_by(Parameter) %>%
    summarise(estimate = mean(value)) %>%
    mutate(method = "Top", calibrated = TRUE)
ens_top_comp <-
    ggs(ens_comp_btl, "p") %>%
    group_by(Parameter) %>%
    summarise(estimate = mean(value)) %>%
    mutate(method = "Compositional", calibrated = TRUE)

csmf_estimates <-
    ens_top_calib %>%
    bind_rows(ens_top_comp)

csmf_estimates <-
    inner_join(csmf_estimates, true_df, by = "Parameter") %>%
    mutate(model = "Ensemble")
csmf_acc_df <-
    csmf_estimates %>%
    group_by(method, calibrated, model) %>%
    summarise(csmfa = csmf_acc(true_value, estimate)) %>%
    mutate(cccsmfa = (csmfa - .632)/(1-.632))

csmf_acc_all <-
    do.call(rbind, posteriors_list) %>%
    bind_rows(csmf_acc_df)

saveRDS(csmf_acc_all, output_file)
quit('no')