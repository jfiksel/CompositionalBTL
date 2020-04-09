library(tidyverse)
library(here)
results <- readRDS(here("simulations", "results_dir.rds"))
csmf_plot <-
    results %>%
    filter(metric == "CSMF") %>%
    ggplot(aes(x = method, y = acc)) +
    geom_boxplot() +
    facet_wrap(~M, nrow = 1) +
    xlab("Method") +
    ylab("CSMF Accuracy")

ccc_plot <-
    results %>%
    filter(metric == "CCC") %>%
    ggplot(aes(x = method, y = acc)) +
    geom_boxplot() +
    facet_wrap(~M, nrow = 1) +
    xlab("Method") +
    ylab("CCC Accuracy")
dir.create(here("simulations", "figs"), recursive = TRUE)
ggsave(here("simulations", "figs", "csmf_acc_dir.pdf"),
       plot = csmf_plot, width = 10, height = 4)
ggsave(here("simulations", "figs", "ccc_acc_dir.pdf"),
       plot = ccc_plot, width = 10, height = 4)