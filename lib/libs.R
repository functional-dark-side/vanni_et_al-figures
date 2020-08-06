nclass.all <- function(x, fun = median)
{
  fun(c(
    nclass.Sturges(x),
    nclass.scott(x),
    nclass.FD(x)
  ))
}

calc_bin_width <- function(x, ...)
{
  rangex <- range(x, na.rm = TRUE)
  (rangex[2] - rangex[1]) / nclass.all(x, ...)
}

get_ratio <- function(x = x, y = y, display = 4/3){
  ratio_display <- display
  ratio_values <- ((max(x) - min(x)))/((max(y) - min(y)))
  ratio_values/ratio_display
}

get_ratio_log10 <- function(x = x, y = y, display = 4/3){
  ratio_display <- display
  ratio_values <- (log10(max(x) - min(x)))/(log10(max(y) - min(y)))
  ratio_values/ratio_display
}

plot_treemap_aa <- function(X){
  library(treemapify)
  X %>%
    mutate(p_pfam_not_covered_aa = pfam_not_covered_aa/n_aa,
           p_pfam_covered_aa = pfam_covered_aa/n_aa) %>%
    select(category_type, cat, pfam_not_covered_aa, pfam_covered_aa ) %>%
    pivot_longer(-contains("cat")) %>%
    ggplot(aes(area = value, subgroup = cat, subgroup2 = category_type, fill = name)) +
    geom_treemap() +
    geom_treemap_subgroup_text(place = "middle", colour = "white", alpha = 0.7, grow = T) +
    geom_treemap_subgroup2_text(colour = "white", alpha = 0.9, fontface = "italic") +
    geom_treemap_subgroup_border(colour = "white", size = 5) +
    geom_treemap_subgroup2_border(colour = "white", size = 2) +
    scale_fill_manual(values = pfam_aa_coverage_colors) +
    theme(legend.position = "none")
}
