library(tidyverse)
library(here)
library(ggpubr)
csmf_acc <- readRDS(here("phmrc_adult_5_causes_analysis", "csmf_phmrc_adult.rds"))
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
    group_by(method, calibrated, model, country, ncalib) %>%
    summarize(csmfa_var = var(csmfa), csmfa = mean(csmfa), cccsmfa = mean(cccsmfa)) %>%
    ungroup()
csmf_acc_mean$method <- recode(csmf_acc_mean$method,
                               Compositional = "Compositional GBQL",
                               Top  = "Single-class/categorical GBQL")
csmf_acc_mean_indiv <-
    csmf_acc_mean %>%
    mutate(ncalib = ifelse(calibrated == "Uncalibrated", 0, ncalib)) %>%
    filter(model != "Ensemble", (calibrated == "Calibrated" | ncalib == 0))

no_labels <- filter(csmf_acc_mean_indiv, ncalib == 0)
no_labels$method <- recode(no_labels$method,
                           `Compositional GBQL` = "PA",
                           `Single-class/categorical GBQL` = "CC")

pa <- filter(no_labels, method == "PA")

acc_apa_df <- readRDS(here("phmrc_adult_5_causes_analysis", "acc_apa.rds"))
acc_apa_df$model <- recode(acc_apa_df$model,
                           insilico = "InSilicoVA",
                           nbc = "NBC",
                           interva = "InterVA",
                           tariff = "Tariff")
apa <- filter(acc_apa_df, method == "APA")


cccsmfa_with_n <-
    csmf_acc_mean_indiv %>%
    mutate(ncalib = ifelse(calibrated == "Uncalibrated", 0, ncalib)) %>%
    filter(model != "Ensemble", (calibrated == "Calibrated" | ncalib == 0)) %>%
    filter(method == "Compositional GBQL") %>%
    ggplot(aes(x = ncalib, y = cccsmfa, color = method)) +
    geom_point() +
    geom_hline(data = pa, aes(yintercept = cccsmfa, color = method)) +
    geom_hline(data = apa, aes(yintercept = cccsmfa, color = method)) +
    geom_line(aes(group = method)) +
    facet_grid(model ~ country) +
    xlab("Number of labeled observations")+
    ylab("CCNAAA") +
    theme(legend.title = element_blank(),
          legend.position = "bottom") +
    scale_color_discrete(breaks=c("PA","APA", "Compositional GBQL"))
figs_dir <- here("phmrc_adult_5_causes_analysis","figs")
if(!dir.exists(figs_dir)){
    dir.create(figs_dir, recursive = TRUE)
}
ggsave(file.path(figs_dir, "cccsmfa_with_n.pdf"),
       cccsmfa_with_n, width = 8, height = 6)

### Now compare compositional GBQL to single cause GBQL
col_scale <- scales::hue_pal()(5)
cols <- c("Ensemble" = col_scale[1],
          "InSilicoVA" = col_scale[2],
          "InterVA" = col_scale[3],
          "NBC" = col_scale[4],
          "Tariff" = col_scale[5])
csmf_acc_single_vs_comp <-
    csmf_acc_mean_indiv %>%
    filter(calibrated == "Calibrated", model != "Ensemble") %>%
    group_by(model, country, ncalib) %>%
    summarize(compositional = cccsmfa[method == "Compositional GBQL"],
              single_cause = cccsmfa[method == "Single-class/categorical GBQL"])
single_vs_comp_plot <-
    ggplot(csmf_acc_single_vs_comp, aes(x = compositional, y = single_cause)) +
    geom_point(aes(color = model)) +
    facet_wrap(~country) +
    scale_color_manual(values = cols) +
    #geom_line(aes(group = method), alpha = .5) +
    facet_wrap( ~ country, nrow = 1) +
    xlab("CCNAA Compositional GBQL")+
    ylab("CCNAA Single-class/categorical \n GBQL") +
    theme(legend.title = element_blank(),
          legend.position = "bottom") +
    geom_abline(slope = 1, alpha = .25)

ggsave(file.path(figs_dir, "cccsmfa_single_vs_comp.pdf"),
       single_vs_comp_plot, width = 8, height = 3)


cccsfma_plot <- ggarrange(cccsmfa_with_n, single_vs_comp_plot,
                          heights = c(2, 1), ncol =1, nrow = 2,
                          labels = c("A", "B"),
                          label.x = .03, label.y = 1)
ggsave(file.path(figs_dir, "cccsmfa_combined.pdf"),
       cccsfma_plot, width = 13, height = 6)





cccsfma_ensemble <-
    csmf_acc_mean %>%
    filter(calibrated == "Calibrated", method != "Top Broad") %>%
    filter(method == "Compositional GBQL") %>%
    ggplot(aes(x = ncalib, y = cccsmfa, color = model)) +
    geom_point() +
    geom_line() +
    facet_wrap( ~ country, nrow = 1) +
    xlab("Number of labeled observations")+
    ylab("CCNAAA") +
    scale_color_manual(values=cols) +
    scale_x_continuous(breaks = seq(0, 400, by = 100),
                       labels = seq(0, 400, by = 100),
                       limits = c(0, 400)) +
    theme(legend.title = element_blank(),
          legend.position = "bottom")
ggsave(file.path(figs_dir, "cccsfma_ensemble.pdf"),
       cccsfma_ensemble, width = 8, height = 3)

# cccsfma_200 <-
#     csmf_acc_mean %>%
#     filter(ncalib == 200, model != "Ensemble", method != "Top Broad") %>%
#     ggplot(aes(x = method, y = cccsmfa)) +
#     geom_point(aes(color = calibrated)) +
#     facet_grid(model ~ country) +
#     xlab("Algorithm Output")+
#     ylab("CCCSMFA") +
#     theme(legend.title = element_blank())
# 
# ggsave(file.path(figs_dir, "cccsmfa_200.pdf"),
#        cccsfma_200, width = 8, height = 8)

# cccsfma_200_ensemble <-
#     csmf_acc_mean %>%
#     filter(ncalib == 200, calibrated == "Calibrated", method != "Top Broad") %>%
#     ggplot(aes(x = model, y = cccsmfa)) +
#     geom_point() +
#     facet_grid(method ~ country) +
#     xlab("Algorithm")+
#     ylab("CCCSMFA") +
#     ylim(.5, .9)+
#     theme(legend.title = element_blank(),
#           axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(file.path(figs_dir, "cccsfma_200_ensemble.pdf"),
#        cccsfma_200_ensemble, width = 6, height = 6)
# 
# 
# countries <- unique(csmf_acc_mean$country)

country_plot_list <- lapply(countries, function(c) {
    p <-
        csmf_acc_mean %>%
        filter(country == c, model != "Ensemble", method != "Top Broad") %>%
        ggplot(aes(x = method, y = cccsmfa)) +
        geom_point(aes(color = calibrated))+
        facet_grid(model ~ ncalib)+
        labs(title = c) +
        xlab("Algorithm Output") +
        ylab("CCCSMFA") +
        theme(legend.title = element_blank())
    ggsave(file.path(figs_dir, paste0("csmf_accuracy_",c,".pdf")),
           p,
           width = 10, height = 4)
    return(p)
})

pdf(file.path(figs_dir, "csmf_accuracy_individual_methods.pdf"), width = 10, height = 4)
for(i in 1:length(country_plot_list)) {
    print(country_plot_list[[i]])
}
dev.off()
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
figs_dir <- here("phmrc_adult_5_causes_analysis","figs")
if(!dir.exists(figs_dir)){
    dir.create(figs_dir, recursive = TRUE)
}
pdf(file.path(figs_dir, "cccsmf_accuracy_individual_methods.pdf"), width = 10, height = 4)
for(i in 1:length(country_plot_list)) {
    print(country_plot_list[[i]])
}
dev.off()

### CCCSMF accuracy

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

### Now same plots for CCC
ccc_acc <- readRDS(here("phmrc_adult_5_causes_analysis", "ccc_phmrc_adult.rds"))
ccc_acc_mean <-
    ccc_acc %>%
    group_by(method, calibrated, model, country, ncalib) %>%
    summarize(ccc = mean(ccc))
countries <- unique(ccc_acc_mean$country)
country_plot_list <- lapply(countries, function(c) {
    p <-
        ccc_acc_mean %>%
        filter(country == c, model != "Ensemble") %>%
        ggplot(aes(x = method, y = ccc)) +
        geom_point(aes(color = calibrated))+
        facet_grid(model ~ ncalib)+
        labs(title = c)
    return(p)
})
figs_dir <- here("phmrc_adult_5_causes_analysis","figs")
if(!dir.exists(figs_dir)){
    dir.create(figs_dir, recursive = TRUE)
}
pdf(file.path(figs_dir, "ccc_accuracy_individual_methods.pdf"), width = 10, height = 4)
for(i in 1:length(country_plot_list)) {
    print(country_plot_list[[i]])
}
dev.off()

### Now make plot to show benefit of ensemble
ccc_acc_ensemble <-
    ccc_acc_mean %>%
    filter(calibrated == TRUE, method == "Compositional") %>%
    ggplot(aes(x = model, y = ccc)) +
    geom_point() +
    facet_grid(country ~ ncalib) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
ggsave(file.path(figs_dir, "ccc_accuracy_ensemble_compositional.pdf"),
       ccc_acc_ensemble, width = 10, height = 4)

