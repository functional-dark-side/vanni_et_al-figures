#!/usr/bin/env Rscript
# Code for Figure 2 panel in Vanni et al.
suppressMessages({
  suppressWarnings({
    # Alluvial plot for the PR results
    library(tidyverse)
    library(ggpubr)
    library(naturalsort)
    library(cowplot)
    library(ggthemr)
    library(RSQLite)

    fig_num <- 2

    db <- "data/Fig2.sqlite"

    con <- RSQLite::dbConnect(RSQLite::SQLite(), db)


    # PANEL C -----------------------------------------------------------------

    cat(paste0("Creating Fig ", fig_num, " - Panel C..."))

    pr_results <- tbl(con, "pr_results") %>% collect()

    pr_results %>%
      group_by(cl_name) %>%
      add_count() %>%
      mutate(cl_name_agg = ifelse(n > 25 || consensus_superkingdom == "Viruses", cl_name, "Other"),
             comb_tax = paste(cl_name_agg, consensus_superkingdom, sep = "#")) %>%
      write_tsv("results/PR_alluvial.tsv")

    cat(" done\n")

    # PANEL D -----------------------------------------------------------------
    cat(paste0("Creating Fig ", fig_num, " - Panel D..."))

    # HQ clusters for the ribosomal proteins
    cl_ribo <- tbl(con, "cl_ribo") %>% collect()
    cl_ribo_hq <- tbl(con, "cl_ribo") %>% collect()

    r_all <- cl_ribo %>%
      select(riboprot, com) %>%
      filter(!is.na(riboprot)) %>%
      unique() %>%
      group_by(riboprot) %>%
      count(name = "n_all")

    r_hq <- cl_ribo_hq %>%
      filter(!is.na(riboprot)) %>%
      select(riboprot, com) %>%
      unique() %>%
      group_by(riboprot) %>%
      count(name = "n_hq")

    ggthemr(layout = "scientific", palette = "fresh")
    fct_sort = function(.f, .fun = sort, ...) {
      f = forcats:::check_factor(.f) # Not needed in dev version
      fct_relevel(f, .fun(levels(f), ...))
    }

    r_comb <- r_all %>%
      ungroup() %>%
      full_join(r_hq) %>%
      replace(is.na(.), 1) %>%
      pivot_longer(-riboprot, names_to = "type", values_to = "counts")

    r_order <- r_comb$riboprot %>% unique() %>% naturalsort(decreasing = TRUE)
    r_plot <- r_comb %>%
      mutate(type = fct_rev(type),
             riboprot = fct_relevel(riboprot, r_order)) %>%
      ggplot(aes(riboprot, counts, fill = type)) +
      geom_col(position = "dodge", color = "black", size = 0.1, width = 0.9) +
      ggpubr::rotate() +
      xlab("Ribosomal proteins") +
      ylab("Number of cluster communities") +
      theme(legend.position = "top") +
      scale_fill_manual(values = c("#E28F4D", "#2A4556"))

    save_plot(filename = "figures/Fig2-hq_cl_riboprot.pdf",  plot = r_plot, base_width = 4, base_height = 4)
    dbDisconnect(con)
    cat(" done\n\nAll figures saved at figures/\n\n")
  })
})
