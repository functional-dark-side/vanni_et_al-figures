library(tidyverse)
library(cowplot)

#Environmental:2,017,658, Genomic:2,347,502, Shared:922,599

data_kept <- tibble(categ = c("Environmental", "Shared", "GTDB r86"),
                    n = c(2017658, 922599, 2347502)) %>%
  mutate(prop = n/sum(n),
         x = "x",
         categ = fct_relevel(categ, rev(c("Environmental", "Shared", "GTDB r86"))))

colors <- c("#D0D5D9", "#CBC6C1", "#F8ECE0")
names(colors) <- c("Environmental", "Shared", "GTDB r86")

p1 <- ggplot(data_kept, aes(x, prop, fill = categ)) +
  geom_col() +
  theme_void() +
  coord_flip() +
  scale_fill_manual(values = colors) +
  theme(aspect.ratio = 0.1,
        legend.position = "top")

save_plot(last_plot(), filename = "figures/Fig1-mg_gtdb_shared_kept.pdf", base_width = 5, base_height = 3)

