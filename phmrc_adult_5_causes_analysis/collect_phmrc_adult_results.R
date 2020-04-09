library(tidyverse)
library(here)
files <- list.files(here("phmrc_adult_5_causes_analysis", "cluster_output"),
                     full.names = TRUE)
seeds <- 1:500
countries <- c("India", "Mexico", "Philippines", "Tanzania")
setting.df <- expand.grid(seed.index = seeds,
                          country = countries,
                          ncalib = c(25, 100, 200, 400))
csmf_df <- do.call(rbind, lapply(files, function(f){
    base <- gsub(".rds", "", basename(f))
    num <- as.numeric(gsub("run-", "", base))
    setting <- setting.df[num,]
    df <- data.frame(readRDS(f), setting)
    return(df)
}))
saveRDS(csmf_df, here("phmrc_adult_5_causes_analysis", "csmf_phmrc_adult.rds"))
#saveRDS(ccc_df, here("phmrc_adult_5_causes_analysis", "ccc_phmrc_adult.rds"))
quit('no')