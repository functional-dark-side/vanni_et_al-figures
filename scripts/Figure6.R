#!/usr/bin/env Rscript
suppressMessages({
  # Phenotype mutant analisys -----------------------------------------------
  library(tidyverse)
  library(gggenes)
  library(ggtree)
  library(ggthemr)
  library(RSQLite)
  library(ape)
  source("lib/colors.R")
  source("lib/libs.R")

  fig_num <- 6

  ggthemr::ggthemr(layout = "scientific", palette = "fresh")

  db <- "data/Fig6.sqlite"
  con <- RSQLite::dbConnect(RSQLite::SQLite(), db)

  # Load data
  # For plotting the Tn-Seq data
  #load("data/fig6/mutants_gu19737823/mutant_genes_cls_gu19737823.Rda", verbose = T)
  mutant_genes_cls <- tbl(con, "mutant_genes_cls_gu19737823") %>% collect()
  sel_genomes <- tbl(con, "sel_genomes_gu19737823") %>% collect()
  sel_genomes_com <- tbl(con, "sel_genomes_com_gu19737823") %>% collect()


  # For plotting the glyph for AO356_08590
  #load("data/fig6/mutants_gu19737823/genes_int_gu19737823.Rda",  verbose = T)
  genes_int <- tbl(con, "genes_int_gu19737823") %>% collect()

  # For plotting the occurrence in OM-RGCv2 and other MG
  #load("data/fig6/mutants_gu19737823/samp_sel_cls_gu19737823.Rda",  verbose = T)
  samp_sel_cls <- tbl(con, "samp_sel_cls_gu19737823") %>% collect()

  #load("data/fig6/mutants_gu19737823/omrgc_genes_gu19737823.Rda",  verbose = T)
  omrgc_genes <- tbl(con, "omrgc_genes_gu19737823") %>% collect()

  gtdb_tax <- read_tsv("data/gtdb/bac_taxonomy_r86.tsv", col_names = c("genome", "taxonomy_string"))
  gtdb_tree <- read.tree("data/gtdb/gtdb_r86_bac.tree")

  # Prepare taxonomy data
  gtdb_tax <- gtdb_tax %>%
    filter(genome %in% gtdb_tree$tip.label) %>%
    mutate(tax_string = gsub("d__|p__|c__|o__|f__|g__|s__", "", taxonomy_string)) %>%
    separate(tax_string,
             into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
             sep = ";",
             remove = TRUE) %>%
    mutate_if(is.character, list(~na_if(.,"")))

  # We load data from Fig3
  db_fig3 <- "data/Fig3.sqlite"

  con_fig3 <- RSQLite::dbConnect(RSQLite::SQLite(), db_fig3)

  # Samples used for the paper
  samples_list <- tbl(con_fig3, "samples_list") %>% collect()

  hmp_cdata <- tbl(con_fig3, "hmp_cdata") %>%
    collect()

  mp_cdata <- tbl(con_fig3, "mp_cdata") %>%
    collect()

  osd_cdata <- tbl(con_fig3, "osd_cdata") %>%
    collect()

  gos_cdata <- tbl(con_fig3, "gos_cdata") %>%
    collect()
  mg_data_filt_by_sample <- tbl(con_fig3, "mg_data_filt_by_sample") %>%
    collect()

  dbDisconnect(con_fig3)

  # Load data for the hhblits graph
  #load("data/fig6/mutants_gu19737823/hhblits_graph_gu_c_12103.Rda", verbose = TRUE)
  hhblits_graph <- tbl(con, "hhblits_graph_gu_c_12103") %>%
    collect() %>%
    .$graph %>%
    unlist(recursive = FALSE) %>%
    unserialize()


  # Load data for the glyphs
  #load("data/fig6/mutants_gu19737823/data_glyphs_gu19737823.Rda", verbose = TRUE)
  data_glyphs_pseudo <- tbl(con, "data_glyphs_pseudo_gu_c_12103") %>% collect()
  data_glyphs_order <- tbl(con, "data_glyphs_order_gu_c_12103") %>% collect()

  # Load data for the tree
  #load("data/fig6/mutants_gu19737823/tree_data_gu19737823.Rda", verbose = TRUE)

  tree_data <- tbl(con, "tree_data_gu19737823") %>% collect()
  gene_tree <- tbl(con, "gene_tree_gu_c_12103") %>%
    collect() %>%
    .$tree %>%
    unlist(recursive = FALSE) %>%
    unserialize()


  # PANEL A & B -------------------------------------------------------------

  cat(paste0("Creating Fig ", fig_num, " - Panel A & B..."))

  # Spectinocym plot for Pseudomonas fluorescens FW300-N2C3
  org <- "pseudo5_N2C3_1"
  base_exp <- "LB"
  tnseq_data <- mutant_genes_cls %>%
    filter(orgId == org) %>% filter(expDesc == base_exp | expDesc == "Spectinomycin 0.025 mg/ml") %>%
    #filter(orgId == org) %>% filter(expDesc == base_exp | (expDesc %in% conds$expDesc)) %>%
    mutate(expDesc = gsub("LB ", "", expDesc)) %>%
    select(locusId, cat, cl_name, expDesc, fit) %>%
    pivot_wider(names_from = expDesc, values_from = fit) %>%
    pivot_longer(cols = c(-contains(base_exp), -contains('locusId'), -contains('cat'), -contains('cl_name')), names_to = "treat", values_to = "fit") %>%
    unite(cl_name, c(cl_name, locusId), sep = ' - ')
  ratio_tnseq_plot <- get_ratio(x = tnseq_data$LB, y = tnseq_data$fit, display = 4/3)

  tnseq_plot <- ggplot(tnseq_data, aes_(as.name(base_exp), ~fit, fill = ~cat, label = ~cl_name)) +
    geom_abline(intercept = 0, size = 0.1, color = "#2F2F2B") +
    geom_hline(yintercept = 0, size = 0.1, color = "#2F2F2B") +
    geom_vline(xintercept = 0, size = 0.1, color = "#2F2F2B") +
    geom_point(shape = 21, alpha = 0.8, color = "#2F2F2B") +
    scale_fill_manual(values = color_comb_cats_I) +
    xlab("Fitness in LB") +
    ylab("Fitness in Spectinomycin 0.025 mg/ml") +
    coord_fixed(ratio = ratio_tnseq_plot) +
    theme(legend.position = "none")


  # Gene glyph for locus: AO356_08590
  plot_gene_int <- ggplot(genes_int, aes(xmin = begin, xmax = end, y = molecule, fill = gene, forward = direction)) +
    geom_gene_arrow() +
    theme_genes() +
    scale_fill_manual(values = gene_colors) +
    theme(legend.position = "top",
          axis.text.x = element_blank(),
          axis.line = element_blank(),
          axis.title = element_blank())


  # Plot occurrence of the interesting cluster in metagenomes
  samples_project <- samples_list %>%
    group_by(project) %>%
    count(name = "samples") %>%
    ungroup()

  samp_sel_cls_plot <- samp_sel_cls %>%
    group_by(project) %>%
    count() %>%
    ungroup() %>%
    inner_join(samples_project) %>%
    mutate(prop = n/samples,
           project = fct_reorder(project, n)) %>%
    ggplot(aes(project, samples)) +
    geom_col(width = 0.7, size = 0.5, color = "#404040", fill = "#545B60", alpha = 0.2) +
    geom_col(aes(project, n), width = 0.7, size = 0.3, color = "black", fill = "#545B60") +
    ggpubr::rotate() +
    scale_y_log10() +
    theme(aspect.ratio = 1/4,
          panel.grid = element_blank()) +
    xlab("") +
    ylab("Metagenomes")


  p1 <- ggpubr::ggarrange(tnseq_plot, samp_sel_cls_plot, grid::nullGrob() )
  cowplot::save_plot(filename = "figures/Fig6-tnseq_plot.pdf",  plot = p1, base_width = 8, base_height = 5)

  cat(" done\n")
  # PANEL C -----------------------------------------------------------------

  cat(paste0("Creating Fig ", fig_num, " - Panel C..."))

  # Graph
  library(tidygraph)
  library(ggraph)
  library(igraph)

  g <- ggraph(hhblits_graph %>%
                mutate(type = ifelse(name == 19737823, "19737823", "Other"),
                       degree = centrality_degree())) +
    geom_edge_fan(aes(color = weight)) +
    geom_node_point(color = "black", aes(fill = type, size = degree), shape = 21) +
    scale_edge_color_distiller(palette = "RdBu", direction = -1) +
    scale_fill_manual(values =  c("#E85A5B", "#37656C"), guide = "none") +
    scale_size(guide = "none")

  write.graph(hhblits_graph %>%
                mutate(type = ifelse(name == 19737823, "19737823", "Other"),
                       degree = centrality_degree()) %>% as.igraph(), file = "results/Fig6-mutant_graph.graphml", format = "graphml")

  cat(" done\n")
  # PANEL D -----------------------------------------------------------------
  cat(paste0("Creating Fig ", fig_num, " - Panel D..."))

  tmp <- gtdb_tax %>%
    as_tibble() %>%
    filter(genome %in% tree_data$label) %>%
    mutate(label = paste0(paste0("o: ", order),"; ", paste0("f: ", family))) %>%
    .$label %>% unique()

  tree_data <- tree_data %>%
    mutate(desc = paste0(paste0(row_number(),"o: ", order),"; ", paste0("f: ", family)))


  cl_counts <- gtdb_tax %>%
    as_tibble() %>%
    filter(genome %in% sel_genomes_com$genome, grepl("GCF", genome)) %>%
    mutate(label = paste0(paste0("o: ", order),"; ", paste0("f: ", family))) %>%
    group_by(order) %>%
    count() %>%
    inner_join(tree_data %>% select(parent, order) %>% group_by(order) %>% arrange(parent) %>% slice(1)) %>%
    ungroup()

  # Plot tree + genes
  all_genes <- bind_rows(data_glyphs_pseudo %>%
                           mutate(gene = case_when(name == "30S ribosomal protein S18" ~ "30S ribosomal protein S18",
                                                   name == "30S ribosomal protein S6" ~ "30S ribosomal protein S6",
                                                   name == "50S ribosomal protein L9" ~ "50S ribosomal protein L9",
                                                   name == "replicative DNA helicase" ~ "replicative DNA helicase",
                                                   gene == "GENE" ~ "gu_c_12103",
                                                   TRUE ~ gene)) %>%
                           mutate(gene = fct_relevel(gene, c("gu_c_12103", "30S ribosomal protein S6", "30S ribosomal protein S18", "50S ribosomal protein L9", "replicative DNA helicase", "OTHER"))),
                         data_glyphs_order %>%
                           mutate(gene = case_when(name == "30S ribosomal protein S18" ~ "30S ribosomal protein S18",
                                                   name == "30S ribosomal protein S6" ~ "30S ribosomal protein S6",
                                                   name == "50S ribosomal protein L9" ~ "50S ribosomal protein L9",
                                                   name == "replicative DNA helicase" ~ "replicative DNA helicase",
                                                   gene == "GENE" ~ "gu_c_12103",
                                                   TRUE ~ gene)) %>%
                           mutate(gene = fct_relevel(gene, c("gu_c_12103", "30S ribosomal protein S6", "30S ribosomal protein S18", "50S ribosomal protein L9", "replicative DNA helicase", "OTHER"))))

  all_genes <- all_genes %>% inner_join(tree_data %>% select(label, desc))

  label2desc <- all_genes %>% select(label, desc) %>% distinct() %>% droplevels()

  gene_tree$tip.label <- plyr::mapvalues(gene_tree$tip.label, from = label2desc$label, to = label2desc$desc)

  to_reverse <- all_genes %>%
    mutate(molecule = desc) %>%
    select(molecule, gene, start, end, strand, direction) %>%
    filter(strand == "reverse", gene == "gu_c_12103") %>% .$molecule
  all_genes_fwd <- all_genes %>%
    mutate(molecule = desc) %>%
    select(molecule, gene, start, end, strand, direction) %>%
    filter(!(molecule %in% to_reverse))
  all_genes_rev <- all_genes %>%
    mutate(molecule = desc) %>%
    select(molecule, gene, start, end, strand, direction) %>%
    filter((molecule %in% to_reverse)) %>%
    group_by(molecule) %>%
    mutate(start = -1 * start, end = -1 * end) %>%
    ungroup()

  p <- ggtree(gene_tree, layout = 'rectangular', aes(color = cl_name)) %<+% (tree_data %>% mutate(label = desc) %>% left_join(cl_counts) %>% mutate(n = ifelse(is.na(n), 0, n))) +
    # geom_tippoint(aes(size = n),
    #               shape = 21,# Make bubbles on edges
    #               fill = "#022641",
    #               color = "#243643",
    #               alpha = 0.7) +
    geom_point2(aes(subset = n > 0, size = n), shape = 21, color = "#454345", fill = "#C0C2C2") +
    geom_tiplab(size = 2.6,
                align = TRUE,
                linesize = 0.2,
                linetype = "dotted",
                color = "black") +
    geom_facet(mapping = aes(xmin = start, xmax = end, fill = gene, forward = direction),
               data = bind_rows(all_genes_fwd, all_genes_rev),
               geom = geom_motif, panel = 'Alignment',
               on = 'gu_c_12103', align = 'left', arrowhead_height = unit(2, "mm"),
               arrowhead_width = unit(1, "mm"), arrow_body_height =  unit(2, "mm")) +
    #scale_color_gradientn(colours = pal, name="Percentage of MAGs", labels=scales::percent) +
    theme(legend.position = "top",
          legend.key = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    scale_fill_manual(values = gene_colors) +
    scale_color_manual(values = cls_colors, na.value = "#454345", guide = "none") +
    scale_size_continuous(range = c(1,6), trans = "sqrt")
  p2 <- facet_widths(p, widths = c(1,2))

  cowplot::save_plot(filename = "figures/Fig6-mutant_tree.pdf",  plot = p2, base_width = 8, base_height = 5)
  dbDisconnect(con)
  cat(" done\nAll figures saved at figures/\n")
})

