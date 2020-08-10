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

StatPercentileX <- ggproto("StatPercentileX", Stat,
                           compute_group = function(data, scales, probs) {
                             percentiles <- quantile(data$x, probs=probs)
                             data.frame(xintercept = percentiles)
                           },
                           required_aes = c("x")
)

stat_percentile_x <- function(mapping = NULL, data = NULL, geom = "vline",
                              position = "identity", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    stat = StatPercentileX, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

StatPercentileXLabels <- ggproto("StatPercentileXLabels", Stat,
                                 compute_group = function(data, scales, probs) {
                                   percentiles <- quantile(data$x, probs=probs)
                                   data.frame(x=percentiles, y=Inf,
                                              label=paste0("p", probs*100, ": ",
                                                           scales::comma(round(10^percentiles, digits=1))))
                                 },
                                 required_aes = c("x")
)

stat_percentile_xlab <- function(mapping = NULL, data = NULL, geom = "text",
                                 position = "identity", na.rm = FALSE,
                                 show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    stat = StatPercentileXLabels, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
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


round_preserve_sum <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

get_stats_aa <- function(X){

  #Total number of aas
  n_aas <- X$n_aa %>% sum()

  aa_summary <- X %>%
    select(cat, contains("_aa"), n_genes) %>%
    group_by(cat) %>%
    summarise(across(everything(), sum))

  k_aa <- aa_summary %>%
    filter(cat == "Known") %>%
    group_by(cat) %>%
    select(cat, contains("_aa")) %>%
    summarise(across(everything(), function(x)x/n_aas))
  unk_aa <- aa_summary %>%
    filter(cat == "Unknown") %>%
    group_by(cat) %>%
    select(cat, contains("_aa")) %>%
    summarise(across(everything(), function(x)x/n_aas))
  nc_aa <- aa_summary %>%
    filter(cat == "None") %>%
    group_by(cat) %>%
    select(cat, contains("_aa")) %>%
    summarise(across(everything(), function(x)x/n_aas))

  stats_prop <- bind_rows(k_aa, unk_aa, nc_aa) %>%
    bind_cols(aa_summary %>% select(n_genes) %>% mutate(n_genes = n_genes/sum(n_genes)))

  stats_prop_wide <- stats_prop %>%
    pivot_longer(-cat) %>%
    mutate(type = ifelse(name == "n_aa", "Total", "Pfam")) %>%
    group_by(type) %>%
    mutate(prop = round_preserve_sum(value, digits = 2)) %>%
    arrange(name)
  list(stats_prop = stats_prop, stats_prop_wide = stats_prop_wide)
}

plot_treemap_aa <- function(X){
  library(treemapify)
  X %>%
    mutate(p_pfam_not_covered_aa = pfam_not_covered_aa/n_aa,
           p_pfam_covered_aa = pfam_covered_aa/n_aa) %>%
    select(category_type, cat, pfam_not_covered_aa, pfam_covered_aa ) %>%
    pivot_longer(-contains("cat")) %>%
    group_by(cat, name) %>%
    summarise(value = sum(value)) %>%
    ggplot(aes(area = value, subgroup = cat, fill = name)) +
    geom_treemap(color = "white", size = 0 , alpha = 1) +
    geom_treemap_subgroup_text(place = "middle", colour = "black", grow = FALSE) +
    geom_treemap_subgroup_border(color = "white", size = 4) +
    scale_fill_manual(values = pfam_aa_coverage_colors) +
    theme_ipsum_rc() +
    theme(legend.position = "none")
}

