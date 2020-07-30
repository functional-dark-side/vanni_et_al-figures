library(tidyverse)
library(tidytree)
library(ape)
library(ggthemr)
library(ggnewscale)
library(cowplot)
library(ggpubr)
library(wesanderson)
library(lvplot)
library(grid)
library(purrr)
library(treemapify)
library(RSQLite)
source("lib/colors.R")
source("lib/libs.R")


db <- "data/Fig5.sqlite"

con <- RSQLite::dbConnect(RSQLite::SQLite(), db)

f1score <- tbl(con, "f1score") %>% collect()

gtdb_phylum_stats_mag_n <- tbl(con, "gtdb_phylum_stats_mag_n")
gtdb_phylum_stats_mag_p <- tbl(con, "gtdb_phylum_stats_mag_p")

gtdb_tree <- read.tree("data/gtdb/gtdb_r86_bac.tree")
gtdb_tax <- read_tsv("data/gtdb/bac_taxonomy_r86.tsv", col_names = c("tip", "taxonomy_string"))

# Prepare taxonomy data
gtdb_tax <- gtdb_tax %>%
  filter(tip %in% gtdb_tree$tip.label) %>%
  mutate(tax_string = gsub("d__|p__|c__|o__|f__|g__|s__", "", taxonomy_string)) %>%
  separate(tax_string,
           into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";",
           remove = TRUE) %>%
  mutate_if(is.character, list(~na_if(.,"")))

# PANEL A -----------------------------------------------------------------
# RED - lineage specific plot

# Load RED data
# Instructions to get RED: https://github.com/functional-dark-side/functional-dark-side.github.io/blob/master/scripts/Phylogenomic_analyses/phylo_plots.R#L207
bac_red_path <- "data/gtdb/gtdb_r86_red/"
bac_red_files <- list.files(path = bac_red_path, pattern = "rank_distribution.tsv", recursive = TRUE, full.names = TRUE)

bac_red_df <- map_dfr(bac_red_files, read_tsv, col_names = T)

bac_red_ranks <- bac_red_df %>% mutate(rank = case_when(grepl("p__", Taxa) ~ "Phylum",
                                                        grepl("c__", Taxa) ~ "Class",
                                                        grepl("o__", Taxa) ~ "Order",
                                                        grepl("f__", Taxa) ~ "Family",
                                                        grepl("g__", Taxa) ~ "Genus",
                                                        grepl("s__", Taxa) ~ "Species")) %>%
  group_by(rank) %>%
  summarise(red = median(`Relative Distance`)) %>%
  ungroup()

ggthemr(layout = "scientific", palette = "fresh")

bac_red_rank_data <- f1score %>%
  group_by(lowest_level, categ) %>%
  count() %>%
  ungroup() %>%
  spread(categ, value = n) %>%
  replace_na(list(EU = 0, GU = 0, KWP = 0, K = 0)) %>%
  mutate(Known = KWP + K, Unknown = EU + GU) %>%
  select(lowest_level, Known, Unknown) %>%
  rename(rank = lowest_level)

bac_red_rank_plot <- bac_red_rank_data %>%
  gather(categ, n, -rank) %>%
  inner_join(bac_red_ranks) %>%
  ggplot(aes(red, n, fill = categ, color = categ, group = categ)) +
  geom_vline(xintercept = bac_red_ranks$red, color = "black", alpha = 0.6, size = 0.3, linetype = "dotdash") +
  #geom_smooth(size=.5) +
  geom_line(size = 0.5) +
  geom_point(shape = 21, color = "black", size = 2, stroke = 0.4, alpha = 0.8) +
  geom_text(aes(x = red, y = 450000, label = rank),
            hjust = 1.1, color = "black") +
  expand_limits(x = 0.25) +
  scale_y_log10(labels = comma) +
  scale_color_manual(values = color_comb_cats_I) +
  scale_fill_manual(values = color_comb_cats_I) +
  ylab("Lineage specific clusters") +
  xlab("Relative evolutionary divergence") +
  coord_flip() +
  theme(panel.grid = element_blank(),
        legend.position = "top")


# PANEL B -----------------------------------------------------------------

##### Phylogenetic conservation (trait depth) ##################################################################################

tD_all <- f1score %>%
  #separate(trait, into = c("cat"), sep = "_", remove = F,extra="drop") %>%
  select(cl_name, cat, mean_depth, var_depth, max_depth, mean_random_depth, P, is_specific) %>%
  mutate(lspec = case_when(is_specific == TRUE ~ "Spec",
                           TRUE ~ "No Spec"))

tD_all_p <- tD_all %>%
  mutate(sign = ifelse(P < 0.05, "sign", "non-sign")) %>%
  group_by(cat) %>%
  add_count() %>%
  rename(tot = n) %>%
  group_by(cat, tot, sign) %>%
  count() %>%
  ungroup() %>%
  mutate(p = n/tot,
         cat = fct_relevel(cat,  c("Unknown", "Known")))

p1 <- ggplot(tD_all_p, aes(cat, p, fill = sign)) +
  geom_col(position = "stack", alpha = 0.9, color = "#404040") +
  #ylab("Proportion of clusters") +
  #xlab("") +
  #geom_boxjitter() +
  scale_y_continuous(labels = percent_format(5L)) +
  #scale_fill_manual(values=c("#E84646","#65ADC2","#556c74", "#233B43")) +
  #scale_color_manual(values=c("black","grey"), labels=c("Non significant (P≥0.05)","Significant (P<0.05)")) +
  scale_fill_manual(values=c("#404040", "#FC3338"), labels=c("Non significant (P ≥ 0.05)","Significant (P < 0.05)")) +
  guides(fill=guide_legend(title = "Non-random distribution")) +
  theme(legend.position = "top") +
  ylab("Proportion of clusters") +
  xlab("")


tD_data <- tD_all %>%
  mutate(prand = ifelse(P < 0.05, "Non-random", "Random")) %>%
  filter(prand == "Non-random")


tD_data_p <- tD_data %>%
  bind_rows(tD_data %>% mutate(lspec = "All"))
stat.test <- compare_means(mean_depth ~ cat, group.by = "lspec", data = tD_data_p)

tau_lv_plot <- ggplot(tD_data_p, aes(cat, mean_depth, fill = cat)) +
  geom_lv(size = 0.5, width.method = "height", color = "#404040", width = 0.5, alpha = 1) +
  stat_compare_means(label = "p.signif", comparisons =  list(c("Known", "Unknown"))) +
  theme(legend.position = "top",
        panel.grid = element_blank()) +
  scale_fill_manual(values = color_comb_cats_I) +
  facet_wrap(~lspec) +
  scale_y_log10() +
  coord_trans(y = 'reverse') +
  ylab(expression("Mean trait depth ("~tau[D]~")")) +
  xlab("")

tD_all %>% mutate(prand = ifelse(P < 0.05, "Non-random", "Random" )) %>%
  group_by(cat, prand) %>%
  count()

tD_all %>% filter(P < 0.05) # non-randomly distributed 465,148

p4 <- ggarrange(grid::nullGrob(), ggarrange(tau_lv_plot, grid::nullGrob(), ncol = 1, nrow = 2), bac_red_rank, ncol = 3, nrow = 1, widths = c(1, 1.5, 2))


save_plot(tau_lv_plot, filename = "figures/Fig5-tau_lv_plot.pdf", base_height = 3.5, base_width = 2.5)

save_plot(plot = p4, filename = "figures/Fig5-gtdb_lineage_mod.pdf", base_height = 5, base_width = 8)
save_plot(plot = quad_plot, filename = "figures/quad_plot.pdf", base_height = 3, base_width = 5)

save_plot(plot = bac_red_rank, filename = "figures/Fig5-bac_red_rank.pdf", base_height = 6.5, base_width = 8.5)



# PANEL C -----------------------------------------------------------------
# Treemap virus -----------------------------------------------------------
f1_cl_phages <- tbl(con, "f1_cl_phages") %>% collect()

phage_data <- f1_cl_phages %>%
  group_by(cl_name, orf_in_phage) %>%
  count() %>%
  ungroup() %>%
  pivot_wider(names_from = orf_in_phage, values_from = n, values_fill = list(n = 0)) %>%
  pivot_longer(-cl_name, names_to = "orf_in_phage", values_to = "n") %>%
  group_by(cl_name) %>%
  mutate(N = sum(n), prop = n/N) %>%
  ungroup() %>%
  filter(orf_in_phage == TRUE, prop > 0.75) %>%
  mutate(has_phage = TRUE)

phage_data <- f1score %>%
  left_join(phage_data %>% select(cl_name, has_phage)) %>%
  mutate(has_phage = ifelse(is.na(has_phage), "No phage", "Phage"),
         is_specific = ifelse(is_specific == TRUE, "Specific", "Non-specific"),
         categ = ifelse(grepl("K", categ), "Known", "Unknown")) %>%
  group_by(is_specific, has_phage, categ) %>%
  count() %>%
  ungroup() %>%
  mutate(has_phage = fct_rev(has_phage))

ggplot(phage_data, aes(area = n, subgroup = is_specific, subgroup2 = has_phage, fill = categ)) +
  geom_treemap() +
  geom_treemap_subgroup_text(place = "middle", colour = "white", alpha = 0.7, grow = T) +
  geom_treemap_subgroup2_text(colour = "white", alpha = 0.9, fontface = "italic") +
  geom_treemap_subgroup_border(colour = "white", size = 2) +
  geom_treemap_subgroup2_border(colour = "white", size = 2) +
  scale_fill_manual(values = color_comb_cats_I) +
  theme(legend.position = "none")

save_plot(plot = last_plot(), filename = "figures/Fig5-phage_treemap.pdf", base_height = 3, base_width = 3)

dbDisconnect(con)


# PANEL E -----------------------------------------------------------------

db <- "data/Fig5-omrgc.sqlite"

con <- RSQLite::dbConnect(RSQLite::SQLite(), db)

cogs_annot <- tbl(con, "omrgc_v2_cogs") %>% collect()
omrgc <- tbl(con, "omrgc_v2_genes_class_categ")

omrgc_gu <- omrgc %>%
  filter(categ == "GU") %>%
  collect()

omrgc_gu_gtdb <- omrgc_gu %>%
  inner_join(f1score %>% select(-categ, -cat))

t_map <- omrgc_gu_gtdb %>%
  mutate(is_specific = ifelse(is_specific == TRUE, TRUE, FALSE),
         lowest_level = as.character(lowest_level),
         lowest_level =  case_when(lowest_level == "Phylum" | lowest_level == "Class" | lowest_level == "Order" ~ "Other",
                                   is.na(lowest_level) ~ "None",
                                   TRUE ~ lowest_level)) %>%
  group_by(is_specific, lowest_level) %>%
  count() %>% ungroup()

ggplot(t_map, aes(area = n, subgroup = is_specific, subgroup2 = lowest_level, fill = lowest_level)) +
  geom_treemap() +
  geom_treemap_subgroup_text(place = "middle", colour = "white", alpha = 0.7, grow = T) +
  geom_treemap_subgroup2_text(colour = "white", alpha = 0.9, fontface = "italic") +
  geom_treemap_subgroup_border(colour = "white", size = 2) +
  geom_treemap_subgroup2_border(colour = "white", size = 2) +
  #scale_fill_manual(values = color_comb_cats_I) +
  theme(legend.position = "none")

save_plot(plot = last_plot(), filename = "figures/Fig5-omrgc_v2_treemap.pdf", base_height = 3, base_width = 3)

dbDisconnect(con)


