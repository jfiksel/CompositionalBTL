library(openVA)
library(caret)
library(here)
library(tidyverse)
library(limSolve)

phmrc <- read.csv(getPHMRC_url("adult"))


### Get COD list
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



country.df <- data.frame(site = c("AP", "Bohol", "Dar", "Mexico", "Pemba", "UP"),
                         country = c("India", 
                                     "Philippines",
                                     "Tanzania",
                                     "Mexico",
                                     "Tanzania",
                                     "India"))
country <- country.df$country[match(phmrc$site, country.df$site)]
models <- c("insilico", "tariff", "interva", "nbc")

setting.df <- expand.grid(country =  unique(as.character(country)), model = models)

### Get CV estimates of M
data.file <- here("phmrc_adult_5_causes_analysis", "M_estimates.rds")
if(!file.exists(data.file)) {
    M_estimates <- lapply(1:nrow(setting.df), function(s) {
        setting <- setting.df[s,]
        c <- setting$country
        model <- setting$model
        train <- phmrc[country != c,]
        #train <- train[1:200,]
        set.seed(123)
        flds <- createFolds(train$gs_text46, k = 10, list = TRUE, returnTrain = TRUE)
        
        M_apa <- M_acc <- matrix(0, nrow = C, ncol = C)
        
        for(k in 1:length(flds)) {
            train.cv <- train[flds[[k]],]
            test.cv <- train[-flds[[k]],]
            ### InSilico
            set.seed(123)
            if(model == 'insilico') {
                phmrc.model <- codeVA(data = test.cv, data.type = "PHMRC",
                                      model = "InSilicoVA",
                                      data.train = train.cv, causes.train = "gs_text34",
                                      phmrc.type = "adult",
                                      jump.scale = 0.05, convert.type = "fixed",
                                      Nsim=10000, auto.length = FALSE)
                ### InSilicoProbs
                insilico_probs <- getIndivProb(phmrc.model)
                insilico_probs[insilico_probs == 0] <- .0000001
                model_probs <- t(apply(insilico_probs, 1, function(x) (x/sum(x))))
            } else if(model == 'tariff') {
                phmrc.model <- codeVA(data = test.cv, data.type = "PHMRC", model = "Tariff",
                                      data.train = train.cv, causes.train = "gs_text34",
                                      phmrc.type = "adult")
                ### Tariff probs
                tariff_score <- phmrc.model$score
                tariff_score[tariff_score == 0] <- .0000001
                ### Take inverse because lower tariff score implies higher degree of belief
                model_probs <- t(apply(tariff_score, 1, function(x) (1/x)/sum(1/x)))
                
            } else if(model == 'interva') {
                phmrc.model <- codeVA(data = test.cv, data.type = "PHMRC", model = "InterVA",
                                      data.train = train.cv, causes.train = "gs_text34",
                                      phmrc.type = "adult")
                ### InterVA Probs
                interva_probs <- getIndivProb(phmrc.model)
                interva_probs[interva_probs == 0] <- .0000001
                model_probs <- t(apply(interva_probs, 1, function(x) (x/sum(x))))
                
            } else {
                interva.probs <- codeVA(data = test.cv, data.type = "PHMRC", model = "InterVA",
                                        data.train = train.cv, causes.train = "gs_text34",
                                        phmrc.type = "adult")
                ### InterVA Probs
                interva_probs <- getIndivProb(interva.probs)
                phmrc.model <- codeVA(data = test.cv, data.type = "PHMRC", model = "NBC",
                                      data.train = train.cv, causes.train = "gs_text34",
                                      phmrc.type = "adult")
                ### NBC probs
                nbc_probs <- getIndivProb(phmrc.model)
                nbc_probs$CaseID <- NULL
                colnames(nbc_probs) <- gsub("\\.", " ", colnames(nbc_probs))
                colnames(nbc_probs)[colnames(nbc_probs) == "Diarrhea Dysentery"] <- "Diarrhea/Dysentery" 
                colnames(nbc_probs)[colnames(nbc_probs) == "Leukemia Lymphomas"] <- "Leukemia/Lymphomas" 
                colnames(nbc_probs)[colnames(nbc_probs) == "Other Non communicable Diseases"] <- "Other Non-communicable Diseases" 
                m <- match(colnames(interva_probs), colnames(nbc_probs))
                nbc_probs <- nbc_probs[,m]
                nbc_probs <- as.matrix(nbc_probs)
                nbc_probs[nbc_probs == 0] <- .0000001
                model_probs <- t(apply(nbc_probs, 1, function(x) (x/sum(x))))
            }
            
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
            true_broad_cause <- cause_map$broad_cause[match(test.cv$gs_text34, cause_map$causes)]
            for(i in 1:C) {
                which.cause <- which(true_broad_cause == causes[i])
                for(j in 1:C) {
                    M_apa[i,j] <- M_apa[i,j] + sum(model_broad_probs[which.cause,j])
                    M_acc[i,j] <- M_acc[i,j] + sum(top_cause_mat[which.cause,j])
                }
            }
        }
        M_apa <- M_apa / rowSums(M_apa)
        M_acc <- M_acc / rowSums(M_acc)
        return(list(M_apa = M_apa, M_acc = M_acc))
    })
    saveRDS(M_estimates, data.file)
} else {
    M_estimates <- readRDS(data.file)
}


### Function to do calibration with set M
calibrationSetM <- function(q, M) {
    ### all arguments should be character vectors
    C <- length(q)
    B <- q
    G <- diag(C)
    H <- rep(0, C)
    E <- matrix(rep(1, C), nrow = 1)
    lsei_fit <- lsei(A = t(M), B = B, G = G, H = H, E = E, F = 1)
    return(lsei_fit$X)
}

### CSMF function
csmf_acc <- function(truth, csmf) {
    acc <- 1 - sum(abs(truth - csmf))/2/(1 - min(truth))
    return(acc)
}

### Collect ACC and APA estimates in list

acc_apa_df <- do.call(rbind, lapply(1:nrow(setting.df), function(s) {
    ### Get setting and M estimates
    setting <- setting.df[s,]
    M_acc <- M_estimates[[s]]$M_acc
    M_apa <- M_estimates[[s]]$M_apa
    
    ### Get CC and PA estimates
    c <- setting$country
    model <- setting$model
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
    cc <- unname(colMeans(top_cause_mat))
    pa <- unname(colMeans(model_broad_probs))
    
    acc_est <- calibrationSetM(cc, M_acc)
    apa_est <- calibrationSetM(pa, M_apa)
    
    test <- phmrc[country == c,]
    true_broad_cause <- cause_map$broad_cause[match(test$gs_text34, cause_map$causes)]
    true_p <- sapply(causes, function(c) mean(true_broad_cause == c))
    
    ### Get CSMF accuracy
    acc_csmfa <- csmf_acc(true_p, acc_est)
    apa_csmfa <- csmf_acc(true_p, apa_est)
    
    df <- data.frame(method = c("ACC", "APA"),
                     calibrated = "Uncalibrated",
                     country = c,
                     ncalib = 0,
                     model = model,
                     csmfa = c(acc_csmfa, apa_csmfa)) %>%
        mutate(cccsmfa = (csmfa - .632)/(1-.632))
    return(df)
}))

### Save data frame
saveRDS(acc_apa_df, here("phmrc_adult_5_causes_analysis", "acc_apa.rds"))

