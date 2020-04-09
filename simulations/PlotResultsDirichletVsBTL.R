library(here)
library(tidyverse)
library(openVA)
results_single_cause <- readRDS(here("simulations", "results_btl_vs_dir_single_cause.rds"))
results_multi_cause <- readRDS(here("simulations", "results_btl_vs_dir_multi_cause.rds"))
results_all <- bind_rows(mutate(results_single_cause, unc = "Known Labels"),
                         mutate(results_multi_cause, unc = "Uncertain Labels"))
results_all$model <- factor(results_all$model,
                                     levels = c("dirichlet-mixture", "overdispered-dirichlet-mixture"),
                                     labels = c("Dirichlet", "Overdispersed Dirichlet"))
results_all$p_U_index <- paste0("p", results_all$p_U_index)
results_all$method <- factor(results_all$method,
                             levels = c("Compositional BTL", "Dirichlet"),
                             labels = c("GBQL", "Dirichlet"))
p_df <- filter(results_all, grepl("p", Parameter))

### First CSMFA plot
csmfa_df <-
    p_df %>%
    group_by(method, model, seed_index, p_U_index, unc) %>%
    summarise(csmfa = getCSMF_accuracy(mean, true_value),
              cccsmfa = (csmfa - .632)/(1-.632)) 

ccnaa_plot <-
    csmfa_df %>%
    filter(unc == "Known Labels") %>%
    group_by(p_U_index, model, method, unc) %>%
    summarize(mean_cccsmfa = mean(cccsmfa)) %>%
    ggplot() +
    geom_point(aes(x = method, y = mean_cccsmfa, color = p_U_index),
               position=position_dodge(width=0.5), size= 3) +
    facet_wrap(~model, ncol = 2) +
    #facet_grid(p_U_index~model) +
    theme(legend.title = element_blank(),
          legend.position = "bottom") +
    ylim(c(0,1)) +
    labs(x = "Model", y = "CCNAA")
figs_dir <- here("simulations", "figs")
ggsave(file.path(figs_dir, "ccnaa_btl_vs_dir.pdf"),
       plot = ccnaa_plot, width = 5, height = 3)


### Rhat and runtime info
mytable <- xtable::xtable(p_df %>%
    filter(unc == "Known Labels") %>%
    filter(model == "Overdispersed Dirichlet") %>%
    group_by(p_U_index, method) %>%
    summarize(avg_rhat = mean(Rhat),
              mean_time = mean(time) / 60) %>%
    pivot_wider(names_from = method, values_from = c(avg_rhat, mean_time)),
    digits =2)
print(mytable, include.rownames = FALSE)

rhat_plot <-
    p_df %>%
    filter(unc == "Known Labels") %>%
    group_by(p_U_index, model, method) %>%
    summarize(mean_rhat = mean(Rhat)) %>%
    ggplot() +
    geom_point(aes(x = method, y = mean_rhat)) +
    facet_grid(p_U_index~model) +
    theme(legend.title = element_blank(),
          legend.position = "bottom") +
    labs(x = "Method", y = "Average Rhat")
figs_dir <- here("simulations", "figs")
ggsave(file.path(figs_dir, "rhat_btl_vs_dir.pdf"),
       plot = rhat_plot, width = 4, height = 4)

runtime_plot <-
    p_df %>%
    filter(unc == "Known Labels") %>%
    group_by(p_U_index, model, method) %>%
    summarize(mean_time = mean(time) / 60) %>%
    ggplot() +
    geom_point(aes(x = method, y = mean_time)) +
    facet_grid(p_U_index~model) +
    theme(legend.title = element_blank(),
          legend.position = "bottom") +
    labs(x = "Method", y = "Average runtime (minutes)")
figs_dir <- here("simulations", "figs")
ggsave(file.path(figs_dir, "time_btl_vs_dir.pdf"),
       plot = runtime_plot, width = 4, height = 4)

#### Compare uncertainty to certainty for GBQL
known_vs_known <-
    csmfa_df %>%
    filter(method == "GBQL") %>%
    group_by(p_U_index, model) %>%
    summarize(ccna_known= mean(cccsmfa[unc == "Known Labels"]),
              ccna_unknown= mean(cccsmfa[unc == "Uncertain Labels"]))
known_vs_unknown_plot <-
    known_vs_known %>%
    ggplot(aes(x = ccna_known, y = ccna_unknown, color = p_U_index, shape=model)) +
    geom_point() +
    geom_abline(slope = 1,col = 'black') +
    theme(legend.title = element_blank(),
          legend.position = "bottom")+
    ylim(.73, .9) +
    xlim(.73, .9) +
    xlab("CCNAA Known Labels") +
    ylab("CCNAA Uncertain Labels ") +
    scale_shape_manual(values = c(17, 15))
ggsave(file.path(figs_dir, "known_vs_unknown_plot.pdf"),
       plot = known_vs_unknown_plot, width = 6, height = 4)


mytable <- xtable::xtable(csmfa_df %>%
    filter(method == "GBQL") %>%
    group_by(p_U_index, model) %>%
    summarize(ccna_known= mean(cccsmfa[unc == "Known Labels"]),
              ccna_unknown= mean(cccsmfa[unc == "Uncertain Labels"])),
    digits = 2)
print(mytable, include.rownames=FALSE)

runtime_boxplot <-
    ggplot(p_single_cause) +
    geom_boxplot(aes(x = method, y = time/60, color = factor(p_U_index))) +
    facet_grid(model ~ Parameter, scales = "free_y") +
    ylab("Time (Minutes)") 



cccsmfa_plot <-
    csmfa_df %>%
    group_by(p_U_index, model, method) %>%
    summarize(mean_cccsmfa = mean(cccsmfa)) %>%
    ggplot() +
    geom_point(aes(x = method, y = mean_cccsmfa))+
    facet_grid(p_U_index~model) +
    labs(y = "CCCSMFA")

m_df <- filter(results_single_cause, grepl("M", Parameter))
cov_m_df <-
    m_df %>%
    group_by(Parameter, method, model) %>%
    summarize(covers = mean(ci_L <= true_value & true_value <= ci_U))

dirichlet_m_plot <-
    m_df %>%
    filter(model == "fitted") %>%
    ggplot() +
    geom_boxplot(aes(x = method, y = mean)) +
    geom_point(aes(x = method, y = true_value), col = 'red') +
    facet_wrap(~ Parameter, nrow = 5) +
    ggtitle("True model dirichlet")

mixture_m_plot <-
    m_df %>%
    filter(model == "mixture-of-mixture") %>%
    mutate(bias = mean - true_value) %>%
    ggplot() +
    geom_boxplot(aes(x = method, y = bias)) +
    geom_point(aes(x = method, y = 0), col = 'red') +
    facet_wrap(~ Parameter, nrow = 5) +
    ggtitle("True model Mixture of Mixtures")

