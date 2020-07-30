#!/usr/bin/env Rscript
# Code for Figure 3 panel in Vanni et al.

library(ggthemr)
library(tidyverse)
library(maditr)
library(RSQLite)
library(cowplot)
library(ggpubr)
library(scales)
source("lib/libs.R")
source("lib/colors.R")
# Panel A -----------------------------------------------------------------

db <- "data/Fig3.sqlite"

con <- RSQLite::dbConnect(RSQLite::SQLite(), db)

# Samples used for the paper
samples_list <- tbl(con, "samples_list") %>% collect()
#filter(study != "OSD", study != "GOS")

# Contextual data
contex <- dbConnect(drv = SQLite(), dbname = "data/contextual_data.db")

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

ggthemr(layout = "scientific", palette = "fresh")
prop_categs_plot <- mg_data_filt_by_sample %>%
  mutate(categ = case_when(categ == "K" ~ "Knowns",
                           categ == "GU" ~ "Genomic unknowns",
                           categ == "EU" ~ "Environmental unknowns",
                           categ == "KWP" ~ "Knowns",
                           TRUE ~ "NC")) %>%
  group_by(categ, biome) %>%
  summarise(abund = sum(abund),
            n_genes = sum(n_genes)) %>%
  group_by(biome) %>%
  mutate(t_n_genes = sum(n_genes),
         t_abund = sum(abund)) %>%
  ungroup() %>%
  mutate(p_n_genes = (n_genes/t_n_genes),
         p_abund = (abund/t_abund)) %>%
  mutate(categ = fct_relevel(categ, color_cats_order_long),
         biome = fct_relevel(biome, c("Marine", "Human"))) %>%
  ggplot(aes(biome, p_n_genes, fill = categ)) +
  geom_col() +
  ggpubr::rotate() +
  scale_fill_manual(values = color_cats_long) +
  scale_y_continuous(labels = percent) +
  ylab("Number of genes") +
  xlab("") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")
save_plot(filename = "figures/Fig3-prop_genes_categ_mg.pdf",  plot = prop_categs_plot, base_width = 5, base_height = 2)
ggthemr_reset()


ggthemr(layout = "scientific", palette = "fresh")
gtdb_props <- tbl(con, "gtdb_props") %>%
  collect() %>%
  mutate(categ = case_when(categ == "KWP" ~ "K",
                           categ == "NONE" ~ "NC",
                           TRUE ~ categ),
         domain = ifelse(domain == "B", "Bacteria", "Archaea")) %>%
  group_by(domain, categ) %>%
  summarise(n = sum(n)) %>%
  group_by(domain) %>%
  mutate(prop = n/sum(n)) %>%
  ungroup()

prop_categs_gtdb_plot <- gtdb_props %>%
  mutate(domain = fct_relevel(domain, c("Bacteria", "Archaea")),
         categ = fct_relevel(categ, color_cats_order)) %>%
  ggplot(aes(domain, prop, fill = categ)) +
  geom_col() +
  ggpubr::rotate() +
  scale_fill_manual(values = color_cats) +
  scale_y_continuous(labels = percent) +
  ylab("Number of genes") +
  xlab("") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")
save_plot(filename = "figures/Fig3-prop_genes_categ_gtdb.pdf",  plot = prop_categs_gtdb_plot, base_width = 5, base_height = 2)
ggthemr_reset()



# Panel B -----------------------------------------------------------------
# Collector curves for metagenomic data. Overall curves
# Without singletons
cum_curve_res_mg_gCl_nosngl <- tbl(con, "cum_curve_res_mg_gCl_nosngl") %>% collect()
# With singletons
cum_curve_res_mg_gCl_sngl <- tbl(con, "cum_curve_res_mg_gCl_sngl") %>% collect()
# Filtered low abundance singletons. Check manuscript for details
cum_curve_res_mg_gCl_sngl_filt <- tbl(con, "cum_curve_res_mg_gCl_sngl_filt") %>% collect()


# Create summaries of the 1000 randomizations
cum_curve_res_mg_gCl_nosngl_summary <- cum_curve_res_mg_gCl_nosngl %>%
  mutate(cat=gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat=c("EU","GU","KWP","K","all"),n=c(0,0,0,0,0),perm=c(0,0,0,0,0),size=c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "Without singletons")

cum_curve_res_mg_gCl_sngl_summary <- cum_curve_res_mg_gCl_sngl %>%
  mutate(cat=gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat=c("EU","GU","KWP","K","all"),n=c(0,0,0,0,0),perm=c(0,0,0,0,0),size=c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "All singletons")

cum_curve_res_mg_gCl_sngl_filt_summary <- cum_curve_res_mg_gCl_sngl_filt %>%
  mutate(cat=gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat=c("EU","GU","KWP","K","all"),n=c(0,0,0,0,0),perm=c(0,0,0,0,0),size=c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "Abundant singletons")

# We used the data with low abundance singletons filtered for Figure 3B
ggthemr(layout = "scientific", palette = "fresh")
cum_curve_plot_gCl_mg <- ggplot() +
  xlab("Metagenomes") +
  ylab("Gene clusters") +
  geom_ribbon(data = cum_curve_res_mg_gCl_sngl_filt_summary %>% filter(cat != "All"), aes(x = size, ymin = mean - sd, ymax = mean + sd, group = interaction(type, cat)), fill = "grey70", alpha = 0.5) +
  geom_line(data = cum_curve_res_mg_gCl_sngl_filt_summary %>% filter(cat != "All") , aes(x = size, y = mean, color = cat, group = interaction(type, cat)), size = 0.7) +
  scale_color_manual(values = color_comb_cats_I) +
  scale_fill_manual(values = color_comb_cats_I) +
  scale_x_continuous(labels = comma, breaks = c(0, 300, 600, 900, 1200)) +
  scale_y_continuous(labels = comma, breaks = c(0, 4e6, 8e6, 12e6)) +
  #scale_linetype_manual(values = c("88", "dashed", "solid")) +
  #scale_x_continuous(limits=c(0,80000000)) +
  theme(legend.position = "top",
        panel.grid = element_blank())

save_plot(filename = "figures/Fig3-col_curve_gCl.pdf",  plot = cum_curve_plot_gCl_mg, base_width = 6, base_height = 6)
ggthemr_reset()

# Supplementart figure comparing all three datasets
cum_curve_res_mg_gCl_summary_comb <- bind_rows(cum_curve_res_mg_gCl_nosngl_summary %>% filter(cat != "All"),
                                               cum_curve_res_mg_gCl_sngl_summary %>% filter(cat != "All"),
                                               cum_curve_res_mg_gCl_sngl_filt_summary %>% filter(cat != "All")) %>%
  mutate(type = fct_relevel(type, c("All singletons", "Abundant singletons", "Without singletons")))


ln_types <- c("dashed", "dotted", "solid")
ggthemr(layout = "scientific", palette = "fresh")
cum_curve_plot_gCl_mg_all <- ggplot() +
  xlab("Metagenomes") +
  ylab("Gene clusters") +
  geom_ribbon(data = cum_curve_res_mg_gCl_summary_comb, aes(x = size, ymin = mean - sd, ymax = mean + sd, group = interaction(type, cat)), fill = "grey70", alpha = 0.5) +
  geom_line(data = cum_curve_res_mg_gCl_summary_comb , aes(x = size, y = mean, color = cat, group = interaction(type, cat), linetype = type), size = 0.7) +
  scale_color_manual(values = color_comb_cats_I) +
  scale_fill_manual(values = color_comb_cats_I) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_linetype_manual(values = ln_types) +
  #scale_x_continuous(limits=c(0,80000000)) +
  theme(legend.position = "top",
        panel.grid = element_blank())

save_plot(filename = "figures/Fig3-col_curve_gCl-sup.pdf",  plot = cum_curve_plot_gCl_mg_all, base_width = 6, base_height = 6)
ggthemr_reset()


# GTDB genomes

# GTDB - genomes ----------------------------------------------------------

cum_curve_res_gtdb_gCl_nosngl <- tbl(con, "cum_curve_res_gtdb_gCl_nosngl") %>% collect()
cum_curve_res_gtdb_gCl_sngl <- tbl(con, "cum_curve_res_gtdb_gCl_sngl") %>% collect()

cum_curve_res_gtdb_gCl_nosngl_summary <- cum_curve_res_gtdb_gCl_nosngl %>%
  mutate(cat=gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat = c("EU","GU","KWP","K","all"),
          n = c(0,0,0,0,0),
          perm = c(0,0,0,0,0),
          size = c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "Without singletons")

cum_curve_res_gtdb_gCl_sngl_summary <- cum_curve_res_gtdb_gCl_sngl %>%
  mutate(cat = gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat = c("EU","GU","KWP","K","all"),
          n = c(0,0,0,0,0),
          perm = c(0,0,0,0,0),
          size = c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "All singletons")

cum_curve_res_gtdb_gCl_summary_comb <- bind_rows(cum_curve_res_gtdb_gCl_nosngl_summary %>% filter(cat != "All"),
                                                 cum_curve_res_gtdb_gCl_sngl_summary %>% filter(cat != "All")) %>%
  mutate(type = fct_relevel(type, c("All singletons", "Without singletons")))


ggthemr(layout = "scientific", palette = "fresh")
ln_types <- c("dashed", "dotted")
cum_curve_plot_gCl_gtdb <- ggplot() +
  xlab("Genomes") +
  ylab("Gene clusters") +
  geom_ribbon(data = cum_curve_res_gtdb_gCl_summary_comb, aes(x = size, ymin = mean - sd, ymax = mean + sd, group = interaction(type, cat)), fill = "grey70", alpha = 0.5) +
  geom_line(data = cum_curve_res_gtdb_gCl_summary_comb , aes(x = size, y = mean, color = cat, group = interaction(type, cat), linetype = type), size = 0.7) +
  scale_color_manual(values = color_comb_cats_I) +
  scale_fill_manual(values = color_comb_cats_I) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_linetype_manual(values = ln_types) +
  theme(legend.position = "top",
        panel.grid = element_blank())

save_plot(filename = "figures/Fig3-col_curve_gCl_gtdb-sup.pdf",  plot = cum_curve_plot_gCl_mg, base_width = 6, base_height = 6)
ggthemr_reset()

ggthemr(layout = "scientific", palette = "fresh")
cum_curve_plot_gCl_gtdb <- ggplot() +
  xlab("Genomes") +
  ylab("Gene clusters") +
  geom_ribbon(data =  cum_curve_res_gtdb_gCl_sngl_summary %>% filter(cat != "All"), aes(x = size, ymin = mean - sd, ymax = mean + sd, group = interaction(type, cat)), fill = "grey70", alpha = 0.5) +
  geom_line(data =  cum_curve_res_gtdb_gCl_sngl_summary %>% filter(cat != "All") , aes(x = size, y = mean, color = cat, group = interaction(type, cat), linetype = type), size = 0.7) +
  scale_color_manual(values = color_comb_cats_I) +
  scale_fill_manual(values = color_comb_cats_I) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  theme(legend.position = "top",
        panel.grid = element_blank())

save_plot(filename = "figures/Fig3-col_curve_gCl_gtdb.pdf",  plot = cum_curve_plot_gCl_gtdb, base_width = 6, base_height = 6)
ggthemr_reset()



# PANEL C -----------------------------------------------------------------
# Metagenomes - biomes ---------------------------------------------------
# MG
cum_curve_res_mg_gCl_filt_HMP <- tbl(con, "cum_curve_res_mg_gCl_filt_HMP") %>% collect()
cum_curve_res_mg_gCl_filt_TARA <- tbl(con, "cum_curve_res_mg_gCl_filt_TARA") %>% collect()
cum_curve_res_mg_gCl_filt_MP <- tbl(con, "cum_curve_res_mg_gCl_filt_MP") %>% collect()
cum_curve_res_mg_gCl_sngl_filt <- tbl(con, "cum_curve_res_mg_gCl_sngl_filt") %>% collect()
cum_curve_res_mg_gCl_filt_marine <- tbl(con, "cum_curve_res_mg_gCl_filt_marine") %>% collect()

cum_curve_res_mg_gCl_filt_HMP_summary <- cum_curve_res_mg_gCl_filt_HMP %>%
  mutate(cat=gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat=c("EU","GU","KWP","K","all"),n=c(0,0,0,0,0),perm=c(0,0,0,0,0),size=c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "HMP")

cum_curve_res_mg_gCl_filt_human_summary <- cum_curve_res_mg_gCl_filt_HMP %>%
  mutate(cat=gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat=c("EU","GU","KWP","K","all"),n=c(0,0,0,0,0),perm=c(0,0,0,0,0),size=c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "Human")

cum_curve_res_mg_gCl_filt_marine_summary <- cum_curve_res_mg_gCl_filt_marine %>%
  mutate(cat=gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat=c("EU","GU","KWP","K","all"),n=c(0,0,0,0,0),perm=c(0,0,0,0,0),size=c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "Marine")

cum_curve_res_mg_gCl_filt_TARA_summary <- cum_curve_res_mg_gCl_filt_TARA %>%
  mutate(cat=gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat=c("EU","GU","KWP","K","all"),n=c(0,0,0,0,0),perm=c(0,0,0,0,0),size=c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "TARA")

cum_curve_res_mg_gCl_filt_MP_summary <- cum_curve_res_mg_gCl_filt_MP %>%
  mutate(cat=gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat=c("EU","GU","KWP","K","all"),n=c(0,0,0,0,0),perm=c(0,0,0,0,0),size=c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "MP")
cum_curve_res_mg_gCl_filt_all_summary <- cum_curve_res_mg_gCl_sngl_filt %>%
  mutate(cat=gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat=c("EU","GU","KWP","K","all"),n=c(0,0,0,0,0),perm=c(0,0,0,0,0),size=c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "All")

cum_curve_res_mg_gCl_summary_comb <- bind_rows(cum_curve_res_mg_gCl_filt_human_summary %>% filter(cat != "All"),
                                               cum_curve_res_mg_gCl_filt_marine_summary %>% filter(cat != "All"))


ln_types <- c("dashed", "dotted", "solid")
ggthemr(layout = "scientific", palette = "fresh")
cum_curve_plot_gCl_mg_biome <- ggplot() +
  xlab("Metagenomes") +
  ylab("Gene clusters") +
  geom_ribbon(data = cum_curve_res_mg_gCl_summary_comb, aes(x = size, ymin = mean - sd, ymax = mean + sd, group = interaction(type, cat)), fill = "grey70", alpha = 0.5) +
  geom_line(data = cum_curve_res_mg_gCl_summary_comb, aes(x = size, y = mean, color = type), size = 0.7) +
  facet_wrap(~cat) +
  scale_color_manual(values = color_biome) +
  scale_fill_manual(values = color_biome) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  #scale_x_continuous(limits=c(0,80000000)) +
  theme(legend.position = "top",
        panel.grid = element_blank())

save_plot(filename = "figures/Fig3-col_curve_gCl-biome.pdf",  plot = cum_curve_plot_gCl_mg_biome, base_width = 6, base_height = 6)
ggthemr_reset()



# SUPPLEMENTARY FIGURES ---------------------------------------------------
# Metagenomes - Viral ---------------------------------------------------
# We also tested the effect of including/exlcuding the TARA viral samples
cum_curve_res_mg_gCl_filt_viral <- tbl(con, "cum_curve_res_mg_gCl_filt_viral") %>% collect()
cum_curve_res_mg_gCl_filt_nonviral <- tbl(con, "cum_curve_res_mg_gCl_filt_nonviral") %>% collect()

cum_curve_res_mg_gCl_filt_viral_summary <- cum_curve_res_mg_gCl_filt_viral %>%
  mutate(cat=gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat=c("EU","GU","KWP","K","all"),n=c(0,0,0,0,0),perm=c(0,0,0,0,0),size=c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "With viral fraction")

cum_curve_res_mg_gCl_filt_nonviral_summary <- cum_curve_res_mg_gCl_filt_nonviral  %>%
  mutate(cat=gsub("-\t","",cat)) %>% ungroup() %>%
  add_row(cat=c("EU","GU","KWP","K","all"),n=c(0,0,0,0,0),perm=c(0,0,0,0,0),size=c(0,0,0,0,0)) %>%
  mutate(cat = case_when(cat == "all" ~ "All",
                         cat == "K" | cat == "KWP" ~ "Known",
                         TRUE ~ "Unknown")) %>%
  group_by(cat, size, perm) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  group_by(cat, size) %>%
  summarise(N = n(),
            mean = mean(n),
            median = median(n),
            min = min(n),
            max = max(n),
            sd = sd(n)) %>%
  mutate(class = "gCl",
         type = "Withot viral fraction")

ggthemr(layout = "scientific", palette = "fresh")
cum_curve_plot_gCl_TARA_nonviral <- ggplot(cum_curve_res_mg_gCl_filt_nonviral_summary %>% filter(cat != "All")) +
  xlab("Metagenomes") +
  ylab("Gene clusters") +
  geom_ribbon(aes(x = size, ymin = mean - sd, ymax = mean + sd, group = interaction(type, cat)), fill = "grey70", alpha = 0.5) +
  #geom_ribbon(data = summary_cum_curve_res_com, aes(x=size,ymin=mean-sd,ymax=mean+sd, group = cat, fill=cat), alpha = 0.3,  color = "grey70", size =.2) +
  #geom_line(data = summary_cum_curve_res_cl, aes(x=size, y=mean, color = cat, group = interaction(class, cat)), size = 0.5) +
  geom_line(aes(x=size, y=mean, color = cat, group = interaction(class, cat)), size = 1) +
  scale_color_manual(values = color_comb_cats_I) +
  scale_fill_manual(values = color_comb_cats_I) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  #scale_x_continuous(limits=c(0,80000000)) +
  theme(legend.position = "top",
        panel.grid = element_blank())
# 10 metagenomes
save_plot(filename = "figures/Fig3-col_curve_gCl_TARA-nonviral-sup.pdf",  plot = cum_curve_plot_gCl_TARA_nonviral, base_width = 6, base_height = 6)
ggthemr_reset()

ggthemr(layout = "scientific", palette = "fresh")
cum_curve_plot_gCl_TARA_viral <- ggplot(cum_curve_res_mg_gCl_filt_viral_summary %>% filter(cat != "All")) +
  xlab("Metagenomes") +
  ylab("Gene clusters") +
  geom_ribbon(aes(x = size, ymin = mean - sd, ymax = mean + sd, group = interaction(type, cat)), fill = "grey70", alpha = 0.5) +
  #geom_ribbon(data = summary_cum_curve_res_com, aes(x=size,ymin=mean-sd,ymax=mean+sd, group = cat, fill=cat), alpha = 0.3,  color = "grey70", size =.2) +
  #geom_line(data = summary_cum_curve_res_cl, aes(x=size, y=mean, color = cat, group = interaction(class, cat)), size = 0.5) +
  geom_line(aes(x=size, y=mean, color = cat, group = interaction(class, cat)), size = 1) +
  scale_color_manual(values = color_comb_cats_I) +
  scale_fill_manual(values = color_comb_cats_I) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  #scale_x_continuous(limits=c(0,80000000)) +
  theme(legend.position = "top",
        panel.grid = element_blank())


save_plot(filename = "figures/Fig3-col_curve_gCl_TARA-viral-sup.pdf",  plot = cum_curve_plot_gCl_TARA_viral, base_width = 6, base_height = 6)
ggthemr_reset()
ggthemr(layout = "scientific", palette = "fresh")
ggarrange(cum_curve_plot_gCl_TARA_viral, cum_curve_plot_gCl_TARA_nonviral, ncol = 2, nrow = 1, common.legend = TRUE, align = "hv")
save_plot(filename = "figures/Fig3-col_curve_gCl_TARA-viral_nonviral-sup.pdf",  plot = last_plot(), base_width = 12, base_height = 6)
ggthemr_reset()


# Calculate the slopes of the curves
get_slopes <- function(X, categ = categ){
  accum <- X %>%
    ungroup() %>%
    filter(cat == categ, size > 0) %>%
    select(cat, size, mean)
  with(accum,diff(mean)/diff(size))
}


get_slopes_avg <- function(X, categ = categ){
  accum <- X %>%
    mutate(cat = gsub("-\t","",cat)) %>%
    ungroup() %>%
    mutate(cat = toupper(cat)) %>%
    mutate(cat = case_when(cat == "ALL" ~ "All",
                           cat == "K" | cat == "KWP" ~ "Known",
                           TRUE ~ "Unknown")) %>%
    group_by(cat, perm, size) %>%
    summarise(n = sum(n)) %>%
    ungroup() %>%
    filter(cat == categ, size > 0)
  perms <- accum$perm %>% unique()
  slopes <- map_dfr(perms, function(Y){
    Z <- accum %>%
      filter(perm == Y) %>%
      arrange(size)
    tibble(slope = with(Z,diff(n)/diff(size)),
           size = Z$size[2:length(Z$size)],
           perm = Y
    )
  })
  slopes %>%
    group_by(size) %>%
    summarise(slope_mean = mean(slope),
              slope_sd = sd(slope))
}

sl_mg_gCl_nosngl_all <- get_slopes(cum_curve_res_mg_gCl_nosngl_summary, "All") %>% enframe(name = "step", value = "diff") %>%  mutate(step = step * 10)
min(sl_mg_gCl_nosngl_all$diff)
sl_mg_gCl_nosngl_k <- get_slopes(cum_curve_res_mg_gCl_nosngl_summary, "Known") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_mg_gCl_nosngl_k$diff)
sl_mg_gCl_nosngl_unk <- get_slopes(cum_curve_res_mg_gCl_nosngl_summary, "Unknown") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_mg_gCl_nosngl_unk$diff)

sl_mg_gCl_sngl_all <- get_slopes(cum_curve_res_mg_gCl_sngl_summary, "All") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_mg_gCl_sngl_all$diff)
sl_mg_gCl_sngl_k <- get_slopes(cum_curve_res_mg_gCl_sngl_summary, "Known") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_mg_gCl_sngl_k$diff)
sl_mg_gCl_sngl_unk <- get_slopes(cum_curve_res_mg_gCl_sngl_summary, "Unknown") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_mg_gCl_sngl_unk$diff)

sl_mg_gCl_sngl_filt_all <- get_slopes(cum_curve_res_mg_gCl_sngl_filt_summary, "All") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_mg_gCl_sngl_filt_all$diff)
sl_mg_gCl_sngl_filt_all %>% filter(diff < 1) %>% arrange(step) %>% head(1)
sl_mg_gCl_sngl_filt_k <- get_slopes(cum_curve_res_mg_gCl_sngl_filt_summary, "Known") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_mg_gCl_sngl_filt_k$diff)
sl_mg_gCl_sngl_filt_k %>% filter(diff < 1) %>% arrange(step) %>% head(1)
sl_mg_gCl_sngl_filt_unk <- get_slopes(cum_curve_res_mg_gCl_sngl_filt_summary, "Unknown") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_mg_gCl_sngl_filt_unk$diff)
sl_mg_gCl_sngl_filt_unk %>% filter(diff < 1) %>% arrange(step) %>% head(1)

# GTDB
sl_gtdb_gCl_nosngl_all <- get_slopes(cum_curve_res_gtdb_gCl_nosngl_summary, "All") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_gtdb_gCl_nosngl_all$diff)
sl_gtdb_gCl_nosngl_k <- get_slopes(cum_curve_res_gtdb_gCl_nosngl_summary, "Known") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_gtdb_gCl_nosngl_k$diff)
sl_gtdb_gCl_nosngl_unk <- get_slopes(cum_curve_res_gtdb_gCl_nosngl_summary, "Unknown") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_gtdb_gCl_nosngl_unk$diff)

sl_gtdb_gCl_sngl_all <- get_slopes(cum_curve_res_gtdb_gCl_sngl_summary, "All") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_gtdb_gCl_sngl_all$diff)
sl_gtdb_gCl_sngl_k <- get_slopes(cum_curve_res_gtdb_gCl_sngl_summary, "Known") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_gtdb_gCl_sngl_k$diff)
sl_gtdb_gCl_sngl_unk <- get_slopes(cum_curve_res_gtdb_gCl_sngl_summary, "Unknown") %>% enframe(name = "step", value = "diff") %>% mutate(step = step * 10)
min(sl_gtdb_gCl_sngl_unk$diff)

p1 <- sl_mg_gCl_nosngl_k %>% mutate(class = "Known", type = "no-sng") %>%
  bind_rows(sl_mg_gCl_nosngl_unk %>% mutate(class = "Unknown", type = "no-sng")) %>%
  bind_rows(sl_mg_gCl_sngl_k %>% mutate(class = "Known", type = "sng")) %>%
  bind_rows(sl_mg_gCl_sngl_unk %>% mutate(class = "Unknown", type = "sng")) %>%
  ggplot(aes(x = step, y = diff, color = class, linetype = type)) +
  geom_line() +
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth() +
  scale_color_manual(values = color_comb_cats_I) +
  scale_fill_manual(values = color_comb_cats_I)


p2 <- sl_gtdb_gCl_nosngl_k %>% mutate(class = "Known", type = "no-sng") %>%
  bind_rows(sl_gtdb_gCl_nosngl_unk %>% mutate(class = "Unknown", type = "no-sng")) %>%
  bind_rows(sl_gtdb_gCl_sngl_k %>% mutate(class = "Known", type = "sng")) %>%
  bind_rows(sl_gtdb_gCl_sngl_unk %>% mutate(class = "Unknown", type = "sng")) %>%
  ggplot(aes(x = step, y = diff, color = class, linetype = type)) +
  #geom_line() +
  scale_y_log10() +
  scale_x_log10() +
  geom_smooth() +
  scale_color_manual(values = color_comb_cats_I) +
  scale_fill_manual(values = color_comb_cats_I)

ggarrange(p1, p2, common.legend = TRUE) + theme_bw()
dbDisconnect(con)
