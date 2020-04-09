library(tidyverse)
library(here)
csmf_acc_known_labels <- readRDS(here("phmrc_adult_5_causes_analysis", "csmf_phmrc_adult.rds"))
csmf_acc_uncertainty <- readRDS(here("phmrc_adult_5_causes_uncertainty", "csmf_phmrc_adult.rds"))
csmf_acc <- bind_rows(mutate(csmf_acc_known_labels, labels = "Known Labels"),
                      mutate(csmf_acc_uncertainty, labels = "Uncertain Labels"))

csmf_acc$model <- factor(csmf_acc$model)
csmf_acc$model <- recode(csmf_acc$model,
                         insilico = "InSilicoVA",
                         nbc = "NBC",
                         interva = "InterVA",
                         tariff = "Tariff")
csmf_acc$calibrated <- ifelse(csmf_acc$calibrated,
                              "Calibrated",
                              "Uncalibrated")
csmf_acc_mean <-
    csmf_acc %>%
    group_by(method, calibrated, model, country, ncalib, labels) %>%
    summarize(csmfa_var = var(csmfa), csmfa = mean(csmfa), cccsmfa = mean(cccsmfa)) %>%
    ungroup() %>%
    filter(method != "Top Broad")
csmf_acc_mean$method <- recode(csmf_acc_mean$method,
                               Compositional = "Compositional GBQL",
                               Top  = "Single-class/categorical GBQL")

csmf_acc_mean <-
    csmf_acc_mean %>%
    mutate(ncalib = ifelse(calibrated == "Uncalibrated", 0, ncalib)) %>%
    filter(model != "Ensemble", (calibrated == "Calibrated" | ncalib == 0)) 

csmf_acc_mean <-
    csmf_acc_mean %>%
    group_by(method, calibrated, model, country, ncalib) %>%
    summarize(n = n(),
              cccsmfa_known = mean(cccsmfa[labels == "Known Labels"]),
              cccsmfa_unknown = mean(cccsmfa[labels == "Uncertain Labels"]))

col_scale <- scales::hue_pal()(5)
cols <- c("InSilicoVA" = col_scale[2],
          "InterVA" = col_scale[3],
          "NBC" = col_scale[4],
          "Tariff" = col_scale[5])
cccsmfa_uncertainty_vs_known <-
    csmf_acc_mean %>%
    filter(method == "Compositional GBQL") %>%
    ggplot(aes(x = cccsmfa_known, y = cccsmfa_unknown,color = model)) +
    geom_point() +
    scale_color_manual(values = cols) +
    #geom_line(aes(group = method), alpha = .5) +
    facet_wrap( ~ country, nrow = 1) +
    xlab("CCNAA Known Labels")+
    ylab("CCNAA Uncertain Labels") +
    theme(legend.title = element_blank(),
          legend.position = "bottom") +
    geom_abline(slope = 1, alpha = .25)
figs_dir <- here("phmrc_adult_5_causes_uncertainty","figs")
if(!dir.exists(figs_dir)){
    dir.create(figs_dir, recursive = TRUE)
}
ggsave(file.path(figs_dir, "cccsmfa_uncertainty_vs_known.pdf"),
       cccsmfa_uncertainty_vs_known, width = 8, height = 3)

### Table
mytable <-
    csmf_acc_mean %>%
    mutate(ncalib = ifelse(calibrated == "Uncalibrated", 0, ncalib)) %>%
    filter(model != "Ensemble", (calibrated == "Calibrated" | ncalib == 0), method == "Probabilistic") %>%

cccsfma_ensemble <-
    csmf_acc_mean %>%
    filter(calibrated == "Calibrated", method != "Top Broad") %>%
    ggplot(aes(x = ncalib, y = cccsmfa, color = model)) +
    geom_point() +
    geom_line() +
    facet_grid(method ~ country) +
    xlab("Number of labeled observations")+
    ylab("CCCSMFA") +
    scale_x_continuous(breaks = seq(0, 400, by = 100),
                       labels = seq(0, 400, by = 100),
                       limits = c(0, 400)) +
    theme(legend.title = element_blank())
ggsave(file.path(figs_dir, "cccsfma_ensemble.pdf"),
       cccsfma_ensemble, width = 10, height = 6)

######


country_plot_list <- lapply(countries, function(c) {
    p <-
        csmf_acc_mean %>%
        filter(country == c, model != "Ensemble") %>%
        ggplot(aes(x = method, y = csmfa)) +
        geom_point(aes(color = calibrated))+
        facet_grid(model ~ ncalib)+
        labs(title = c) +
        ylim(.5, 1)
    return(p)
})
figs_dir <- here("phmrc_adult_5_causes_uncertainty","figs")
if(!dir.exists(figs_dir)){
    dir.create(figs_dir, recursive = TRUE)
}
pdf(file.path(figs_dir, "csmf_accuracy_individual_methods.pdf"), width = 10, height = 4)
for(i in 1:length(country_plot_list)) {
    print(country_plot_list[[i]])
}
dev.off()

### CCCSMF accuracy
country_plot_list <- lapply(countries, function(c) {
    p <-
        csmf_acc_mean %>%
        filter(country == c, model != "Ensemble") %>%
        ggplot(aes(x = method, y = cccsmfa)) +
        geom_point(aes(color = calibrated))+
        facet_grid(model ~ ncalib)+
        labs(title = c) +
        ylim(-.15, 1)
    return(p)
})
figs_dir <- here("phmrc_adult_5_causes_uncertainty","figs")
if(!dir.exists(figs_dir)){
    dir.create(figs_dir, recursive = TRUE)
}
pdf(file.path(figs_dir, "cccsmf_accuracy_individual_methods.pdf"), width = 10, height = 4)
for(i in 1:length(country_plot_list)) {
    print(country_plot_list[[i]])
}
dev.off()



### Look at variance
country_plot_list <- lapply(countries, function(c) {
    p <-
        csmf_acc_mean %>%
        filter(country == c, model != "Ensemble", calibrated==TRUE) %>%
        ggplot(aes(x = method, y = csmfa_var)) +
        geom_point()+
        facet_grid(model ~ ncalib)+
        labs(title = c)
    return(p)
})
figs_dir <- here("phmrc_adult_5_causes_analysis","figs")
if(!dir.exists(figs_dir)){
    dir.create(figs_dir, recursive = TRUE)
}
pdf(file.path(figs_dir, "variance_csmf_accuracy_individual_methods.pdf"), width = 10, height = 4)
for(i in 1:length(country_plot_list)) {
    print(country_plot_list[[i]])
}
dev.off()

### Now make plot to show benefit of ensemble
csmf_acc_ensemble <-
    csmf_acc_mean %>%
    filter(calibrated == TRUE, method == "Compositional") %>%
    ggplot(aes(x = model, y = csmfa)) +
    geom_point() +
    facet_grid(country ~ ncalib) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave(file.path(figs_dir, "csmf_accuracy_ensemble_compositional.pdf"),
       csmf_acc_ensemble, width = 10, height = 4)

