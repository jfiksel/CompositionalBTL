output_dir1 <- "cluster_output_btl_vs_dir_single_cause"
output_dir2 <- "cluster_output_btl_vs_dir_multi_cause"

files1 <- list.files(output_dir1, full.names = TRUE)
results_df1 <- do.call(rbind, lapply(files1, readRDS))
saveRDS(results_df1, "results_btl_vs_dir_single_cause.rds")

files2 <- list.files(output_dir2, full.names = TRUE)
results_df2 <- do.call(rbind, lapply(files2, readRDS))
saveRDS(results_df2, "results_btl_vs_dir_multi_cause.rds")
quit('no')