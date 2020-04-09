library(openVA)
if(!dir.exists("../data")){
    dir.create("../data", recursive = TRUE)
}
if(!file.exists("../data/phmrc_adult.rds")){
    phmrc <- read.csv(getPHMRC_url("adult"))
    saveRDS(phmrc, "../data/phmrc_adult.rds") 
} else {
    phmrc <- readRDS("../data/phmrc_adult.rds")
}
if(!dir.exists('adult_models')){
    dir.create('adult_models', recursive = TRUE)
}

country.df <- data.frame(site = c("AP", "Bohol", "Dar", "Mexico", "Pemba", "UP"),
                         country = c("India", 
                                     "Philippines",
                                     "Tanzania",
                                     "Mexico",
                                     "Tanzania",
                                     "India"))
country <- country.df$country[match(phmrc$site, country.df$site)]

### Train all countries
for(c in unique(as.character(country))) {
    test <- phmrc[country == c,]
    train <- phmrc[country != c,]
    if(!file.exists(paste0("adult_models/insilico_model_", c, ".rds"))){
        set.seed(123)
        phmrc.insilicova <- codeVA(data = test, data.type = "PHMRC",
                                   model = "InSilicoVA",
                                   data.train = train, causes.train = "gs_text34",
                                   phmrc.type = "adult",
                                   jump.scale = 0.05, convert.type = "fixed",
                                   Nsim=10000, auto.length = FALSE)
        saveRDS(phmrc.insilicova, paste0("adult_models/insilico_model_", c, ".rds"))
    } else {
        phmrc.insilicova <- readRDS(paste0("adult_models/insilico_model_", c, ".rds"))
    }
    if(!file.exists(paste0("adult_models/tariff_model_", c, ".rds"))){
        set.seed(123)
        tariff <- codeVA(data = test, data.type = "PHMRC", model = "Tariff",
                         data.train = train, causes.train = "gs_text34",
                         phmrc.type = "adult")
        saveRDS(tariff, paste0("adult_models/tariff_model_", c, ".rds"))
    } else {
        tariff <- readRDS(paste0("adult_models/tariff_model_", c, ".rds"))
    } 
    if(!file.exists(paste0("adult_models/interva_model_", c, ".rds"))){
        set.seed(123)
        interva <- codeVA(data = test, data.type = "PHMRC", model = "InterVA",
                         data.train = train, causes.train = "gs_text34",
                         phmrc.type = "adult")
        saveRDS(interva, paste0("adult_models/interva_model_", c, ".rds"))
    } else {
        interva <- readRDS(paste0("adult_models/interva_model_", c, ".rds"))
    } 
    if(!file.exists(paste0("adult_models/nbc_model_", c, ".rds"))){
        set.seed(123)
        nbc <- codeVA(data = test, data.type = "PHMRC", model = "NBC",
                          data.train = train, causes.train = "gs_text34",
                          phmrc.type = "adult")
        saveRDS(nbc, paste0("adult_models/nbc_model_", c, ".rds"))
    } else {
        nbc <- readRDS(paste0("adult_models/nbc_model_", c, ".rds"))
    } 
    
    ### InSilicoProbs
    insilico_probs <- getIndivProb(phmrc.insilicova)
    insilico_probs[insilico_probs == 0] <- .0000001
    insilico_probs <- t(apply(insilico_probs, 1, function(x) (x/sum(x))))
    saveRDS(insilico_probs, paste0("adult_models/insilico_model_", c, "_probs.rds"))
    
    ### Tariff probs
    tariff_score <- tariff$score
    tariff_score[tariff_score == 0] <- .0000001
    ### Take inverse because lower tariff score implies higher degree of belief
    tariff_probs <- t(apply(tariff_score, 1, function(x) (1/x)/sum(1/x)))
    saveRDS(tariff_probs, paste0("adult_models/tariff_model_", c, "_probs.rds"))
    
    ### InterVA Probs
    interva_probs <- getIndivProb(interva)
    interva_probs[interva_probs == 0] <- .0000001
    interva_probs <- t(apply(interva_probs, 1, function(x) (x/sum(x))))
    saveRDS(interva_probs, paste0("adult_models/interva_model_", c, "_probs.rds"))
    
    ### NBC probs
    nbc_probs <- getIndivProb(nbc)
    nbc_probs$CaseID <- NULL
    colnames(nbc_probs) <- gsub("\\.", " ", colnames(nbc_probs))
    colnames(nbc_probs)[colnames(nbc_probs) == "Diarrhea Dysentery"] <- "Diarrhea/Dysentery" 
    colnames(nbc_probs)[colnames(nbc_probs) == "Leukemia Lymphomas"] <- "Leukemia/Lymphomas" 
    colnames(nbc_probs)[colnames(nbc_probs) == "Other Non communicable Diseases"] <- "Other Non-communicable Diseases" 
    m <- match(colnames(interva_probs), colnames(nbc_probs))
    nbc_probs <- nbc_probs[,m]
    nbc_probs <- as.matrix(nbc_probs)
    nbc_probs[nbc_probs == 0] <- .0000001
    nbc_probs <- t(apply(nbc_probs, 1, function(x) (x/sum(x))))
    saveRDS(nbc_probs, paste0("adult_models/nbc_model_", c, "_probs.rds"))
}
