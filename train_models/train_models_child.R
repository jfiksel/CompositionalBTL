library(openVA)
if(!dir.exists("../data")){
    dir.create("../data", recursive = TRUE)
}
if(!file.exists("../data/phmrc_child.rds")){
    ### Read in child data
    child.raw <- read.csv(getPHMRC_url("child"))
    saveRDS(child.raw, "../data/phmrc_child.rds") 
} else {
    child.raw <- readRDS("../data/phmrc_child.rds")
}
if(!dir.exists('child_models')){
    dir.create('child_models', recursive = TRUE)
}
### Get data frame matching va34 code to actual COD
set.seed(123)
cause.df <- unique(child.raw[,c("gs_text34", "va34")])
cause.df$va34 <- as.character(as.numeric(cause.df$va34))
### Clean data into output usable by Tariff & InsilicoVA
child.clean <- ConvertData.phmrc(child.raw, phmrc.type = "child")$output
### Assign countries
countries <- rep(NA, nrow(child.raw))
countries[child.raw$site %in% c("AP", "UP")] <- "India"
countries[child.raw$site %in% c("Mexico")] <- "Mexico"
countries[child.raw$site %in% c("Dar", "Pemba")] <- "Tanzania"
countries[child.raw$site %in% c("Bohol")] <- "Philippines"
child.clean$Cause <- as.character(cause.df$gs_text34[match(child.clean$Cause, cause.df$va34)])
saveRDS(child.clean, "../data/phmrc_child_clean.rds")
saveRDS(countries, "../data/phmrc_child_countries.rds")

### Train all countries
for(c in c('India', 'Tanzania')) {
    country.data <- child.clean[countries == c,]
    train.final <- child.clean[countries != c,]
    if(!file.exists(paste0("child_models/insilico_model_", c, ".rds"))){
        set.seed(123)
        phmrc.insilicova <- codeVA(data = country.data,
                                   data.type = "customize", model = "InSilicoVA",
                                   data.train = train.final, causes.train = "Cause",
                                   jump.scale = 0.05, Nsim=10000,
                                   auto.length = FALSE)
        saveRDS(phmrc.insilicova, paste0("child_models/insilico_model_", c, ".rds"))
    } else {
        phmrc.insilicova <- readRDS(paste0("child_models/insilico_model_", c, ".rds"))
    }
    if(!file.exists(paste0("child_models/tariff_model_", c, ".rds"))){
        set.seed(123)
        tariff <- codeVA(data = country.data,
                         data.type = "customize", model = "Tariff",
                         data.train = train.final, causes.train = "Cause")
        saveRDS(tariff, paste0("child_models/tariff_model_", c, ".rds"))
    } else {
        tariff <- readRDS(paste0("child_models/tariff_model_", c, ".rds"))
    } 
    if(!file.exists(paste0("child_models/interva_model_", c, ".rds"))){
        set.seed(123)
        interva <- codeVA(data = country.data, data.type = "customize",
                          model = "InterVA",
                          data.train = train.final, causes.train = "Cause")
        saveRDS(interva, paste0("child_models/interva_model_", c, ".rds"))
    } else {
        interva <- readRDS(paste0("child_models/interva_model_", c, ".rds"))
    } 
    if(!file.exists(paste0("child_models/nbc_model_", c, ".rds"))){
        set.seed(123)
        nbc <- codeVA(data = country.data, data.type = "customize", model = "NBC",
                      data.train = train.final, causes.train = "Cause")
        saveRDS(nbc, paste0("child_models/nbc_model_", c, ".rds"))
    } else {
        nbc <- readRDS(paste0("child_models/nbc_model_", c, ".rds"))
    } 
    
    ### InSilicoProbs
    insilico_probs <- getIndivProb(phmrc.insilicova)
    insilico_probs[insilico_probs == 0] <- .0000001
    insilico_probs <- t(apply(insilico_probs, 1, function(x) (x/sum(x))))
    saveRDS(insilico_probs, paste0("child_models/insilico_model_", c, "_probs.rds"))
    
    ### Tariff probs
    tariff_score <- tariff$score
    tariff_score[tariff_score == 0] <- .0000001
    ### Take inverse because lower tariff score implies higher degree of belief
    tariff_probs <- t(apply(tariff_score, 1, function(x) (1/x)/sum(1/x)))
    saveRDS(tariff_probs, paste0("child_models/tariff_model_", c, "_probs.rds"))
    
    ### InterVA Probs
    interva_probs <- getIndivProb(interva)
    interva_probs[interva_probs == 0] <- .0000001
    interva_probs <- t(apply(interva_probs, 1, function(x) (x/sum(x))))
    saveRDS(interva_probs, paste0("child_models/interva_model_", c, "_probs.rds"))
    
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
    saveRDS(nbc_probs, paste0("child_models/nbc_model_", c, "_probs.rds"))
}
