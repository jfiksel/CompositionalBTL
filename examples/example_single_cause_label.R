library(here)
library(gtools)
library(here)
library(ggmcmc)
library(coda)
library(tidyverse)
library(stringr)
library(MCMCpack)
### If you want to use here, this should work
### May have to initiate .Rproj for CompositionalBTL
source(here("GibbsSamplingScripts", "CompositionalBTLGibbs.R"))
source(here("GibbsSamplingScripts", "CompositionalBTLGibbs_pshrink.R"))
### Otherwise you can set the working directory to current working directory and run
# source(file.path("..", "GibbsSamplingScripts", "CompositionalBTLGibbs.R"))

### Will subset to test data = India, ncalib = 200
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
set.seed(123)
seeds <- sample(-1e6:1e6, size = 500, replace = F)
init.seed = seeds[seed.index]
set.seed(init.seed)
### Sample labeled individuals using inverse probabilities
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

### Show posterior draws for InSilicoVA
k <- 1
insilico_comp_btl <- compositional_btl(A_U = comp_array[,,k],
                                       A_L = comp_array[calib_indices,,k],
                                       G_L = G_L,
                                       causes = as.character(1:C),
                                       thin = 5,
                                       burnin = 1000, ndraws = ndraw_gbql,
                                       power = 1/100,
                                       epsilon = .01,
                                       alpha = 10)

### Get calibrated predictions
insilico_calib_p <-
    ggs(insilico_comp_btl, "p") %>%
    group_by(Parameter) %>%
    summarise(estimate = mean(value)) 

insilico_comp_btl_uncalib <- compositional_btl(A_U = comp_array[,,k],
                                       causes = as.character(1:C),
                                       thin = 5,
                                       burnin = 1000, ndraws = ndraw_gbql,
                                       power = 1/100,
                                       epsilon = .01,
                                       alpha = 10)

insilico_calib_p_uncalib <-
    ggs(insilico_comp_btl2, "p") %>%
    group_by(Parameter) %>%
    summarise(estimate = mean(value)) 

colMeans(comp_array[,,k]) - as.vector(insilico_calib_p_uncalib$estimate)

### calibration using p-shrinkage
lambdavec=10^(seq(-3,2,length=50))
insilico_calib_pshrink_mat=Reduce('rbind',lapply(lambdavec,function(lambda){
print(lambda)
insilico_comp_btl_pshrink <- compositional_btl_pshrink(A_U = comp_array[,,k],
                                       A_L = comp_array[calib_indices,,k],
                                       G_L = G_L,
                                       lambda=lambda,
                                       causes = as.character(1:C),
                                       thin = 5,
                                       burnin = 1000, ndraws = ndraw_gbql,
                                       power = 1/100,
                                       epsilon = .01)

param_df <- ggs(insilico_comp_btl_pshrink)
p_df <- ggs(insilico_comp_btl_pshrink, family = "p")
rhat_max <- max(ggs_Rhat(param_df)$data$Rhat)
rhat_p_max = max(ggs_Rhat(p_df)$data$Rhat)

insilico_calib_pshrink <-
    ggs(insilico_comp_btl_pshrink, "p") %>%
    group_by(Parameter) %>%
    summarise(estimate = mean(value)) %>%
    as.data.frame() %>%
    mutate(lambda=lambda,Rhat=rhat_max,Rhat_p=rhat_p_max)
    insilico_calib_pshrink
}))

true_p_df=data.frame(Parameter=paste0("p[",1:5,"]"),estimate=true_p[causes])
true_q_df=data.frame(Parameter=paste0("p[",1:5,"]"),estimate=colMeans(comp_array[,,k]))

### shrinkage paths
insilico_calib_pshrink_mat %>%
    ggplot(aes(x=log10(lambda),y=estimate,col=Parameter)) +
    geom_line() +
    geom_point(data=true_p_df,aes(x=2,y=estimate,col=Parameter),alpha=0.5,size=2) +
    geom_point(data=true_q_df,aes(x=2,y=estimate,col=Parameter),shape=2,size=2)
    
csmfvec=sapply(lambdavec,function(l) 
    csmf_acc(true_p,(insilico_calib_pshrink_mat %>% filter(lambda==l))$estimate))
dev.new()
plot(log10(lambdavec),csmfvec)

rhatvec=sapply(lambdavec,function(l) unique((insilico_calib_pshrink_mat %>% filter(lambda==l))$Rhat))
dev.new()
plot(log10(lambdavec),pmin(rhatvec,5))
abline(h=c(1.01,1.05,1.1))

rhatpvec=sapply(lambdavec,function(l) unique((insilico_calib_pshrink_mat %>% filter(lambda==l))$Rhat_p))
dev.new()
plot(log10(lambdavec),pmin(rhatpvec,5))
abline(h=c(1.01,1.05,1.1))

### both shrinkage paths for estimates and csmf show sharp zigzagging changes but stabilizes after a while
### running for lambda with best csmf
insilico_comp_btl_pshrink_best <- compositional_btl_pshrink(A_U = comp_array[,,k],
                                       A_L = comp_array[calib_indices,,k],
                                       G_L = G_L,
                                       lambda=lambdavec[which.max(csmfvec)],
                                       causes = as.character(1:C),
                                       thin = 5,
                                       burnin = 1000, ndraws = ndraw_gbql,
                                       power = 1/100,
                                       epsilon = .01)

param_df <- ggs(insilico_comp_btl_pshrink_best)
p_df <- ggs(insilico_comp_btl_pshrink_best, family = "p")
rhat_max <- max(ggs_Rhat(param_df)$data$Rhat)
rhat_p_max = max(ggs_Rhat(p_df)$data$Rhat)



### Now ensemble
ens_comp_btl <- compositional_ensemble_btl(A_U = comp_array,
                                           A_L = comp_array[calib_indices,,],
                                           G_L = G_L,
                                           causes = as.character(1:C),
                                           thin = 5,
                                           burnin = 1000, ndraws = ndraw_gbql,
                                           power = 1/100,
                                           epsilon = .01,
                                           alpha = 10)
ens_calib_p <-
    ggs(ens_comp_btl, "p") %>%
    group_by(Parameter) %>%
    summarise(estimate = mean(value)) %>%
    mutate(method = "Compositional", calibrated = TRUE)

csmf_acc(true_p,ens_calib_p$estimate)
q_ens=as.vector(apply(comp_array,2,mean))
csmf_acc(true_p,q_ens)

### calibration using p-shrinkage
lambdavecens=10^(seq(-2,3,length=50))
ensemble_calib_pshrink_mat=Reduce('rbind',lapply(lambdavecens,function(lambda){
print(lambda)
ensemble_comp_btl_pshrink <- compositional_ensemble_btl_pshrink(A_U = comp_array,
                                       A_L = comp_array[calib_indices,,],
                                       G_L = G_L,
                                       lambda=lambda,
                                       causes = as.character(1:C),
                                       thin = 5,
                                       burnin = 1000, ndraws = ndraw_gbql,
                                       power = 1/100,
                                       epsilon = .01)

param_df <- ggs(ensemble_comp_btl_pshrink)
p_df <- ggs(ensemble_comp_btl_pshrink, family = "p")
rhat_max <- max(ggs_Rhat(param_df)$data$Rhat)
rhat_p_max = max(ggs_Rhat(p_df)$data$Rhat)

ensemble_calib_pshrink <-
    ggs(ensemble_comp_btl_pshrink, "p") %>%
    group_by(Parameter) %>%
    summarise(estimate = mean(value)) %>%
    as.data.frame() %>%
    mutate(lambda=lambda,Rhat=rhat_max,Rhat_p=rhat_p_max)
    ensemble_calib_pshrink
}))

true_q_ens_df=data.frame(Parameter=paste0("p[",1:5,"]"),estimate=q_ens)

### shrinkage paths
ensemble_calib_pshrink_mat %>%
    ggplot(aes(x=log10(lambda),y=estimate,col=Parameter)) +
    geom_line() +
    geom_point(data=true_p_df,aes(x=3,y=estimate,col=Parameter),alpha=0.5,size=2) +
    geom_point(data=true_q_ens_df,aes(x=3,y=estimate,col=Parameter),shape=2,size=2)
    
csmfvec=sapply(lambdavecens,function(l) 
    csmf_acc(true_p,(ensemble_calib_pshrink_mat %>% filter(lambda==l))$estimate))
dev.new()
plot(log10(lambdavecens),csmfvec)

rhatvec=sapply(lambdavecens,function(l) unique((ensemble_calib_pshrink_mat %>% filter(lambda==l))$Rhat))
dev.new()
plot(log10(lambdavecens),pmin(rhatvec,5))
abline(h=c(1.01,1.05,1.1))

rhatpvec=sapply(lambdavecens,function(l) unique((ensemble_calib_pshrink_mat %>% filter(lambda==l))$Rhat_p))
dev.new()
plot(log10(lambdavecens),pmin(rhatpvec,5))
abline(h=c(1.01,1.05,1.1))