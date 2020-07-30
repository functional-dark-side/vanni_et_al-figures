#!/usr/bin/env Rscript
library(tidyverse)
library(maditr)
library(RSQLite)
library(cowplot)
library(ggpubr)
library(grid)
library(ggthemr)
library(scales)
source("lib/libs.R")
source("lib/colors.R")
# Figure

# We load data from Fig3
db <- "data/Fig3.sqlite"

con <- RSQLite::dbConnect(RSQLite::SQLite(), db)

# Samples used for the paper
samples_list <- tbl(con, "samples_list") %>% collect()

hmp_cdata <- tbl(con, "hmp_cdata") %>%
  collect()

mp_cdata <- tbl(con, "mp_cdata") %>%
  collect()

osd_cdata <- tbl(con, "osd_cdata") %>%
  collect()

gos_cdata <- tbl(con, "gos_cdata") %>%
  collect()


mg_data_filt_by_sample <- tbl(con, "mg_data_filt_by_sample") %>%
  collect()

data <- mg_data_filt_by_sample %>%
  filter(!(label %in% gos_cdata$label)) %>%
  let(class = case_when(grepl("^K", categ) ~ "Known",
                        TRUE ~ "Unknown"),
      categ = case_when(grepl("^K", categ) ~ "K",
                        TRUE ~ categ)) %>%
  take(abund = sum(abund),
       n_genes = sum(n_genes),
       by = c("type", "categ", "label", "biome", "sample_abund", "sample_size", "class")) %>%
  let(p_abund = abund/sample_abund, p_n_genes = n_genes/sample_size) %>%
  as_tibble() %>%
  filter(type != "DISC") %>%
  mutate(categ = case_when(categ == "K" ~ "Known",
                           categ == "GU" ~ "Genomic unknown",
                           categ == "EU" ~ "Environmental unknown"),
         categ = fct_relevel(categ, c("Known", "Genomic unknown", "Environmental unknown")),
         type = case_when(type == "SINGL" ~ "Singletons",
                          type == "LT10" ~ "< 10 genes",
                          type == "G10" ~ ">= 10 genes"),
         type = fct_relevel(type, c("Singletons", "< 10 genes", ">= 10 genes")))

dbDisconnect(con)

db <- "data/Fig4.sqlite"

con <- RSQLite::dbConnect(RSQLite::SQLite(), db)

# PANEL A -----------------------------------------------------------------
ratio_cat_plot <- get_ratio(x = data$p_n_genes, y = data$p_abund)

cat_plot <- ggplot(data, aes(p_n_genes, p_abund, fill  = biome, label = label)) +
  geom_point(shape = 21, alpha = 0.8, color = "#404040", size = 2) +
  scale_fill_manual(values = color_biome) +
  scale_color_manual(values = color_biome) +
  scale_x_continuous(labels = percent, limits = c(0, 0.82)) +
  scale_y_continuous(labels = percent) +
  facet_grid(categ~type) +
  #coord_fixed(ratio = ratio_cat_plot) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "none") +
  xlab("Gene count") +
  ylab("Gene abundance")
save_plot(filename = "figures/Fig4-categ_props.pdf",  plot = cat_plot, base_width = 10, base_height = 8)


# PANEL B -----------------------------------------------------------------

# HMP
hmp_data <- mg_data_filt_by_sample %>%
  let(class = case_when(grepl("^K", categ) ~ "Known",
                        TRUE ~ "Unknown"),
      categ = case_when(grepl("^K", categ) ~ "K",
                        TRUE ~ categ)) %>%
  take(abund = sum(abund),
       n_genes = sum(n_genes),
       by = c("label", "sample_abund", "sample_size", "categ")) %>%
  let(p_abund = abund/sample_abund, p_n_genes = n_genes/sample_size) %>%
  as_tibble() %>%
  select(label, categ, p_abund) %>%
  filter(categ != "NONE") %>%
  pivot_wider(names_from = categ, values_from = p_abund) %>%
  mutate(rGU = GU/K, rEU = EU/K, rK = K/(EU+GU)) %>%
  inner_join(hmp_cdata) %>%
  mutate(HMP_isolation_body_site = case_when(HMP_isolation_body_site == "airways" ~ "Airways",
                                             HMP_isolation_body_site == "gastrointestinal_tract" ~ "Gastrointestinal tract",
                                             HMP_isolation_body_site == "oral" ~ "Oral",
                                             HMP_isolation_body_site == "skin" ~ "Skin",
                                             HMP_isolation_body_site == "urogenital_tract" ~ "Urogenital tract"))

ratio_hmp <-  get_ratio(x = hmp_data$rEU, y = hmp_data$rGU)
hmp_plot <- ggplot(hmp_data, aes(rEU, rGU, fill = HMP_isolation_body_site, label = label)) +
  #stat_density2d(color = "black", alpha = 0.6) +
  geom_point(shape = 21,  color = "black", alpha = 0.8, size = 1.5) +
  theme_bw() +
  scale_fill_manual(values = color_body_site, name = NULL) +
  #coord_fixed(ratio = ratio_hmp) +
  scale_x_log10() +
  scale_y_log10() +
  theme(panel.grid = element_blank()) +
  ylab("[Genomic unknowns] / [Knowns]") +
  xlab("[Environmental unknowns] / [Knowns]")

save_plot(filename = "figures/Fig4-hmp_plot.pdf",  plot = hmp_plot, base_width = 6, base_height = 4)

crassphage_data <- tbl(con, "crassphage_data") %>% collect()
hpv_data <- tbl(con, "hpv_data") %>% collect()

hmp_plot_crassphage <- ggplot(hmp_data, aes(rEU, rGU, fill = HMP_isolation_body_site, label = label)) +
  #stat_density2d(color = "black", alpha = 0.6) +
  geom_point(shape = 21, fill = "grey",  color = "black", alpha = 0.4, size = 1.5) +
  geom_point(data = hmp_data %>% filter(label %in% crassphage_data$label), shape = 21,  color = "black", alpha = 0.9, size = 1.5) +
  theme_bw() +
  scale_fill_manual(values = color_body_site, name = NULL) +
  #coord_fixed(ratio = ratio_hmp) +
  scale_x_log10() +
  scale_y_log10() +
  theme(panel.grid = element_blank()) +
  ylab("[Genomic unknowns] / [Knowns]") +
  xlab("[Environmental unknowns] / [Knowns]") + ggtitle("Crassphages")

hmp_plot_hpv <- ggplot(hmp_data, aes(rEU, rGU, fill = HMP_isolation_body_site, label = label)) +
  #stat_density2d(color = "black", alpha = 0.6) +
  geom_point(shape = 21, fill = "grey",  color = "black", alpha = 0.4, size = 1.5) +
  geom_point(data = hmp_data %>% filter(label %in% hpv_data$label), shape = 21,  color = "black", alpha = 0.9, size = 1.5) +
  theme_bw() +
  scale_fill_manual(values = color_body_site, name = NULL) +
  #coord_fixed(ratio = ratio_hmp) +
  scale_x_log10() +
  scale_y_log10() +
  theme(panel.grid = element_blank()) +
  ylab("[Genomic unknowns] / [Knowns]") +
  xlab("[Environmental unknowns] / [Knowns]") + ggtitle("HPV")

ggarrange(hmp_plot_crassphage, hmp_plot_hpv, common.legend = TRUE)
save_plot(filename = "figures/Fig4_crassphage_hpv_plot-sup.pdf",  plot = last_plot(), base_width = 10, base_height = 6)



# PANEL C -----------------------------------------------------------------
tara_data <- mg_data_filt_by_sample %>%
  let(class = case_when(grepl("^K", categ) ~ "Known",
                        TRUE ~ "Unknown"),
      categ = case_when(grepl("^K", categ) ~ "K",
                        TRUE ~ categ)) %>%
  take(abund = sum(abund),
       n_genes = sum(n_genes),
       by = c("label", "sample_abund", "sample_size", "categ")) %>%
  let(p_abund = abund/sample_abund, p_n_genes = n_genes/sample_size) %>%
  as_tibble() %>%
  select(label, categ, p_abund) %>%
  filter(categ != "NONE") %>%
  pivot_wider(names_from = categ, values_from = p_abund) %>%
  mutate(rGU = GU/K, rEU = EU/K, rK = K/(EU+GU)) %>%
  filter(grepl("TARA", label)) %>%
  mutate(filter_fraction = case_when(grepl("<-0.22", label)    ~ "Virus-enriched",
                                     grepl("0.1-0.22", label)  ~ "Girus/A|B-enriched",
                                     grepl("0.22-0.45", label) ~ "Girus/A|B-enriched",
                                     grepl("0.45-0.8", label)  ~ "A|B-enriched",
                                     grepl("0.22-1.6", label)  ~ "A|B-enriched",
                                     grepl("0.22-3", label)    ~ "A|B-enriched")) %>%
  mutate(depth_category = case_when( grepl("SRF", label)  ~ 'SRF',
                                     grepl("DCM", label) ~ 'DCM',
                                     grepl("MES", label) ~ 'MES',
                                     grepl("MIX", label) ~ 'MIX'))

ratio_tara <-  get_ratio(x = tara_data$rEU, y = tara_data$rGU, display = 0.5)
tara_plot <- ggplot(tara_data, aes(rEU, rGU, fill = filter_fraction)) +
  #stat_density2d(color = "black", alpha = 0.6) +
  geom_point(shape = 21,  color = "black", alpha = 0.8, size = 1.5) +
  theme_bw() +
  #coord_fixed(ratio = ratio_tara) +
  scale_x_log10() +
  scale_y_log10() +
  theme(panel.grid = element_blank()) +
  ylab("[Genomic unknowns] / [Knowns]") +
  xlab("[Environmental unknowns] / [Knowns]")  +
  scale_fill_manual(name = "", values = color_tara_filter)
save_plot(filename = "figures/Fig4-tara_plot.pdf",  plot = tara_plot, base_width = 6, base_height = 4)

p1 <- plot_grid(hmp_plot + theme(legend.position = "none", text = element_text(size = 11)), tara_plot + theme(legend.position = "none", text = element_text(size = 11)), ncol = 1, nrow = 2, align = "hv")

p2 <- plot_grid(cat_plot + theme(legend.position = "none", text = element_text(size = 11)), plot_grid(nullGrob(), p1, ncol = 1, nrow = 2, rel_heights = c(0.1,0.9)), nrow = 1, ncol = 2, rel_widths = c(1.8, 1), rel_heights = c(2, 1))

save_plot(filename = "figures/Fig4-hmp_tara_plot.pdf",  plot = p2, base_width = 10, base_height = 6)


# PANEL D -----------------------------------------------------------------


# Niche breadth -----------------------------------------------------------
cl_nb_all_mv <- tbl(con, "cl_nb_all_mv") %>%
  collect()

cat_order <- c("Knowns", "Genomic unknowns", "Environmental unknowns")
sign_order <- c("Narrow", "Non significant", "Broad")
cl_nb_all_mv_summary <- cl_nb_all_mv %>%
  mutate(categ_c = case_when(grepl("K", categ) ~ "Knowns",
                             categ == "GU" ~ "Genomic unknowns",
                             TRUE ~ "Environmental unknowns")) %>%
  select(categ_c, sign_mv) %>%
  group_by(categ_c, sign_mv) %>%
  count() %>%
  group_by(categ_c) %>%
  mutate(prop = n/sum(n)) %>%
  ungroup() %>%
  mutate(categ_c = fct_relevel(categ_c, rev(cat_order)),
         sign_mv = fct_relevel(sign_mv, (sign_order)))

dust <- ggthemr(layout = "scientific", palette = "dust", set_theme = FALSE)
ggthemr(layout = "scientific", palette = "greyscale", set_theme = TRUE)
p_nb_cl <- ggplot(cl_nb_all_mv_summary, aes(categ_c, prop, fill = sign_mv)) +
  geom_col(width = 0.7, color = "#404040") +
  ggpubr::rotate() +
  scale_y_continuous(labels = percent) +
  dust$scales$scale_fill_discrete(name = "") +
  xlab("") +
  ylab("Proportion") +
  theme(legend.position = "top",
        aspect.ratio = 1/7,
        legend.key.size = unit(3,"mm"))



com_nb_all_mv <- tbl(con, "com_nb_all_mv") %>% collect()

cat_order <- c("Knowns", "Genomic unknowns", "Environmental unknowns")
sign_order <- c("Narrow", "Non significant", "Broad")
com_nb_all_mv_summary <- com_nb_all_mv %>%
  mutate(categ_c = case_when(grepl("K", categ) ~ "Knowns",
                             categ == "GU" ~ "Genomic unknowns",
                             TRUE ~ "Environmental unknowns")) %>%
  select(categ_c, sign_mv) %>%
  group_by(categ_c, sign_mv) %>%
  count() %>%
  group_by(categ_c) %>%
  mutate(prop = n/sum(n)) %>%
  ungroup() %>%
  mutate(categ_c = fct_relevel(categ_c, rev(cat_order)),
         sign_mv = fct_relevel(sign_mv, (sign_order)))

p_nb_com <- ggplot(com_nb_all_mv_summary, aes(categ_c, prop, fill = sign_mv)) +
  geom_col(width = 0.7, color = "#404040") +
  ggpubr::rotate() +
  scale_y_continuous(labels = percent) +
  dust$scales$scale_fill_discrete(name = "") +
  xlab("") +
  ylab("Proportion") +
  theme(legend.position = "top",
        aspect.ratio = 1/7,
        legend.key.size = unit(3,"mm"))

p3 <- ggarrange(p_nb_cl, p_nb_com, nrow = 2, ncol = 1, common.legend = TRUE)

p4 <- ggarrange(p2, p3, nrow = 2, heights = c(1, 0.5))

save_plot(filename = "figures/Fig4-panelsABCD.pdf",  plot = p4, base_width = 14, base_height = 6)

dbDisconnect(con)


