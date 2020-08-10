library(tidyverse)
library(RSQLite)

source("lib/libs.R")


# Initialize DB
db <- "data/filter-mg-samples.sqlite"
con <- RSQLite::dbConnect(RSQLite::SQLite(), db)

# Get assemblies
HMP1_I_assm <- tbl(con, "hmp1_I_assm") %>%
  collect()

# Read gene data
gene_data <- tbl(con, "gene_data") %>%
  collect()

# Read HMP QC samples
HMP_qc <- tbl(con, "hmp_qc") %>%
  collect()

# Get all cdata
# HMP1-I date = 2011
# HMP1-II date = 2014
HMP_cdata <- tbl(con, "hmp_cdata") %>%
  collect()


# Samples that had a problem during mapping to the assemblies
HMP_nocoverage <- tbl(con, "hmp_nocoverage") %>%
  collect()

# Get bad HMP1-I samples
# 745 samples
HMP1_I <- HMP_cdata %>%
  inner_join(HMP1_I_assm)

# 690 HQ
HMP_bad <- HMP1_I %>%
  filter(phase == "HMP1-I") %>%
  filter(!(label %in% HMP_qc$label))


# Filter out bad samples
gene_data_filt <- gene_data %>%
  filter(!(label %in% HMP_bad$label))

gene_data %>%
  group_by(study) %>%
  count()

#filtering out bad HMP
gene_data_filt %>%
  group_by(study) %>%
  count()

nsamples <- gene_data_filt %>%
  group_by(study) %>%
  count()

gene_data_filt_summary <- summary(gene_data_filt$n_genes)

ggthemr::ggthemr(palette = "fresh", layout = "scientific")

gene_data_filt %>%
  ggplot(aes(n_genes)) +
  geom_histogram(binwidth = calc_bin_width(log10(gene_data_filt$n_genes)), color = "black", alpha = 0.7) +
  stat_percentile_x(probs = c(0.25, 0.5, 0.75), linetype = 2, color = "#B7144B") +
  stat_percentile_xlab(probs = c(0.25, 0.5, 0.75), hjust = 1, vjust = 1.5, angle = 90, size = 2.2) +
  scale_x_log10(labels = scales::comma) +
  xlab("Number of genes") +
  ylab("Counts") +
  theme(axis.title = element_text(size=7),
        axis.text = element_text(size=6))

# We will keep samples with >= genes than the ones in the first quartile
p25 <- quantile(gene_data_filt$n_genes, probs=0.25)

geneXsample <- gene_data_filt %>%
  filter(n_genes >= p25) %>%
  group_by(study) %>%
  count(sort = TRUE)

# Get final set of samples !!!! ##
samples_list <- gene_data_filt %>%
  filter(n_genes >= p25) %>%
  filter(!(label %in% HMP_nocoverage$label))

dbDisconnect(con)
