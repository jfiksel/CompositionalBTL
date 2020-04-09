output_dir <- "cluster_output_dir"
files <- list.files(output_dir, full.names = TRUE)
results_df <- do.call(rbind, lapply(files, readRDS))
saveRDS(results_df, "results_dir.rds")
quit('no')