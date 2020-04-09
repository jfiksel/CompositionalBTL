library(tidyverse)
library(here)
files <- list.files(here("phmrc_adult_5_causes_uncertainty", "cluster_output"),
                     full.names = TRUE)
csmf_df <- lapply(files, readRDS)
csmf_df <- do.call(rbind, csmf_df)
saveRDS(csmf_df, here("phmrc_adult_5_causes_uncertainty", "csmf_phmrc_adult.rds"))
quit('no')
