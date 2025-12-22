############################################################
## Genome-wide coverage profiles for MPXV-S05 across pools
## Pools: RUN-3, RUN-6, RUN-9
## Requires:
##   results-raw-3/MPXV-S05.depth.tsv
##   results-raw-6/MPXV-S05.depth.tsv
##   results-raw-9/MPXV-S05.depth.tsv
############################################################

## Load libraries
library(readr)
library(dplyr)
library(purrr)
library(ggplot2)

##----------------------------------------------------------
## 1. User settings
##----------------------------------------------------------

## Sample ID (used in filenames and plot title)
#sample_id <- "MPXV-S05"
sample_id<-"MPXV-S05"
sample_id2 <-"MPXV-S06"

## Directories for each pooling level
run_dirs1   <- c("results-raw-3", "results-raw-6", "results-raw-9")
run_labels1 <- c("RUN-3",         "RUN-6",         "RUN-9")

## Directories for each pooling level
run_dirs2   <- c("results-raw-6", "results-raw-9")
run_labels2 <- c("RUN-6",         "RUN-9")

## Output PDF name
out_pdf <- paste0("coverage_", sample_id, "_RUN9_RUN15.pdf")

## Bin size in bp
bin_size <- 10000

## set working directory manually if needed
setwd("/path/to/VSP-Training/Manuscript/BMJ-Submission/Revision-01")

##----------------------------------------------------------
## 2. Helper: read one depth file and tag with run label
##----------------------------------------------------------
## Assumes samtools depth format: chrom  pos  depth

read_depth_one <- function(run_dir, run_label, sample_id) {
  file_path <- file.path(run_dir, paste0(sample_id, ".depth.tsv"))
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  read_tsv(
    file_path,
    col_names = c("chrom", "pos", "depth"),
    show_col_types = FALSE
  ) %>%
    mutate(run = run_label)
}

##----------------------------------------------------------
## 3. Read all runs into one data frame
##----------------------------------------------------------

depth_all <- map2_dfr(run_dirs, run_labels, read_depth_one, sample_id = sample_id)

## Quick sanity check
message("Loaded depth data:")
print(dplyr::glimpse(depth_all))

##----------------------------------------------------------
## 4. Bin coverage in 10 kb windows and compute mean depth
##----------------------------------------------------------

depth_binned <- depth_all %>%
  mutate(
    bin_start = (pos %/% bin_size) * bin_size,
    bin_mid   = bin_start + bin_size / 2
  ) %>%
  group_by(run, chrom, bin_start, bin_mid) %>%
  summarise(mean_depth = mean(depth, na.rm = TRUE), .groups = "drop") %>%
  mutate(run = factor(run, levels = c("RUN-3", "RUN-6", "RUN-9")))


depth_binned <- depth_all %>%
  mutate(
    bin_start = (pos %/% bin_size) * bin_size,
    bin_mid   = bin_start + bin_size / 2
  ) %>%
  group_by(run, chrom, bin_start, bin_mid) %>%
  summarise(mean_depth = mean(depth, na.rm = TRUE), .groups = "drop") %>%
  mutate(run = factor(run, levels = c("RUN-6", "RUN-9")))

##----------------------------------------------------------
## 5. Plot genome-wide coverage profiles
##----------------------------------------------------------

## Adjust colours to match your existing BMJ figures if needed
run_cols <- c(
  "RUN-3" = "#1b9e77",  # green
  "RUN-6" = "#d95f02",  # orange
  "RUN-9" = "#7570b3"   # purple
)

p1 <- ggplot(depth_binned,
            aes(x = bin_mid/1000, y = mean_depth, colour = run)) +
  geom_line(alpha = 0.6, linewidth = 0.6) +
  geom_smooth(se = FALSE, linewidth = 0.9, span = 0.15) +
  scale_colour_manual(values = run_cols, name = "Pooling level") +
  labs(
    title = paste("Genome-wide coverage for", sample_id, "across pooling levels"),
    x = "Genomic position (kbp, 10 kb bins)",
    y = "Mean depth"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "top"
  )

## Show plot 
print(p1)


p1+p2

library(cowplot)

# composite plot
plot_grid(p1, p2,
          ncol = 1,       # Arrange in one column (stacked vertically)
          align = "v",    # Specify vertical alignment
          axis = "l"      # Optional: Align by the left y-axis (useful for different axis labels)
)
##----------------------------------------------------------
## 6. Save PDF
##----------------------------------------------------------

ggsave(
  filename = out_pdf,
  plot     = p,
  width    = 7,
  height   = 4,
  units    = "in"
)

message("Saved PDF: ", out_pdf)






############################################################
## MPXV Taxonomic Classification – Final Script
## Author: Julien A. Nguinkal
## Date: 12-10-2025
## Purpose:
##   1. Load Excel table with MPXV classification %
##   2. Clean column names (janitor)
##   3. Remove NA rows
##   4. Convert to long format
##   5. Plot two figures:
##        - Option B  : Grouped bar plot (mean +/- SD)
##        - Option D  : Distribution plot (box + jitter)
##   6. Save high-quality PDFs
############################################################

##----------------------------------------------------------
## 0. Load required libraries
##----------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(janitor)

##----------------------------------------------------------
## 1. Read & clean dataset
##----------------------------------------------------------

# File path 
input_file <- "Taxonomicclassification.xlsx"

tax_raw <- read_excel(input_file) %>%
  clean_names()

cat("\n== Raw data after clean_names ==\n")
print(tax_raw)

# EXPECTED COLUMNS:
# pools, kraken_uniq_percent, kraken2_percent, centrifuge_percent

##----------------------------------------------------------
## 2. Remove rows where all classifiers are NA
##----------------------------------------------------------

tax_clean <- tax_raw %>%
  filter(!(is.na(kraken_uniq_percent) &
             is.na(kraken2_percent) &
             is.na(centrifuge_percent)))

cat("\n== After removing NA rows ==\n")
print(tax_clean)

## Ensure pool levels are factors in correct biological order
tax_clean <- tax_clean %>%
  mutate(
    pools = factor(pools, levels = c(3, 6, 9, 15))
  )

##----------------------------------------------------------
## 3. Reshape to long format & convert to numeric
##----------------------------------------------------------

tax_long <- tax_clean %>%
  pivot_longer(
    cols = c(kraken_uniq_percent, kraken2_percent, centrifuge_percent),
    names_to = "classifier",
    values_to = "percent_mpxv"
  ) %>%
  mutate(
    classifier = recode(
      classifier,
      "kraken_uniq_percent" = "KrakenUniq",
      "kraken2_percent"     = "Kraken2",
      "centrifuge_percent"  = "Centrifuge"
    ),
    classifier   = factor(classifier, levels = c("Kraken2","KrakenUniq","Centrifuge")),
    percent_mpxv = as.numeric(percent_mpxv)
  ) %>%
  filter(!is.na(percent_mpxv))  # Remove remaining NA values

cat("\n== Long format ==\n")
print(head(tax_long, 10))

##----------------------------------------------------------
## 4. Summary table for bar plots: mean ± SD
##----------------------------------------------------------

tax_summary <- tax_long %>%
  group_by(pools, classifier) %>%
  summarise(
    mean_percent = mean(percent_mpxv, na.rm = TRUE),
    sd_percent   = sd(percent_mpxv, na.rm = TRUE),
    n            = n(),
    .groups = "drop"
  )

cat("\n== Summary table (mean ± SD) ==\n")
print(tax_summary)

##----------------------------------------------------------
## 5. Colour palette
##----------------------------------------------------------

class_cols <- c(
  "Kraken2"     = "#1b9e77",
  "KrakenUniq"  = "#d95f02",
  "Centrifuge"  = "#7570b3"
)


##----------------------------------------------------------
## 6. Option Distribution plot (box + jitter)
##----------------------------------------------------------

p_dist <- ggplot(tax_long,
                 aes(x = pools, y = percent_mpxv, fill = classifier)) +
  geom_boxplot(
    position = position_dodge(width = 0.7),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    aes(colour = classifier),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
    size = 1.3,
    alpha = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = class_cols, name = "Classifier") +
  scale_colour_manual(values = class_cols) +
  scale_y_continuous(limits = c(50, 100)) +   # zoom from 50 to 100
  labs(
    title = "",
    x = "Pooling level (samples per run)",
    y = "MPXV-classified reads (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0)
  )


p_dist +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 14),
    axis.text  = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0)
  )


print(p_dist)

ggsave("Figure_tax_classification_distribution.pdf",
       p_dist, width = 7, height = 4, units = "in")

############################################################
## END OF SCRIPT
############################################################



###############################################################################
# FULL VISUAL QC PIPELINE FOR MANUSCRIPT FIGURES
# Author: Julien A. Nguinkal
# Goal: Produce clean plots for QC and variant metrics,
#       grouped purely by Pool Size (3, 6, 9, 15).
###############################################################################

library(tidyverse)
library(janitor)
library(ggpubr)

set.seed(123)   # reproducibility of jitter when few points

###############################################################################
# 1. LOAD INPUT DATA
###############################################################################

# Variant metrics
vars_raw <- read.delim("~/VSP-Training/Manuscript/BMJ-Sumission/variants_data.csv")

# QC metrics
qc_raw   <- read.delim("~/VSP-Training/Manuscript/BMJ-Sumission/data_qc.csv")


###############################################################################
# 2. CLEAN VARIANT DATA
###############################################################################

vars <- vars_raw %>%
  clean_names() %>%                               # pool_id, total_variants, ts_tv, ...
  filter(!is.na(pool_id)) %>%                     # keep only rows with pool assignment
  mutate(
    Pool = factor(pool_id, levels = sort(unique(pool_id))),
    TotalVariants = total_variants,
    TsTv = ts_tv
  ) %>%
  select(Pool, TotalVariants, TsTv)

print(names(vars))
###############################################################################
# 3. CLEAN QC COVERAGE DATA
###############################################################################

qc <- qc_raw %>%
  clean_names() %>%                                # pool_size, genome_ge_10x_pct, mean_depth…
  filter(!is.na(pool_size)) %>%                    # remove empty rows
  mutate(
    Pool      = factor(pool_size, levels = sort(unique(pool_size))),
    Pct10x    = genome_ge_10x_pct,
    MeanDepth = mean_depth
  ) %>%
  select(Pool, Pct10x, MeanDepth)

print(names(qc))


###############################################################################
# 4. FIGURE A — DISTRIBUTION OF VARIANT COUNT BY POOL SIZE
###############################################################################

pA <- vars %>%
  ggplot(aes(x = Pool, y = TotalVariants, fill = Pool)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, color = "black", alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic(base_size = 14) +
  labs(
    title = "Total Variant Count by Pool Size",
    x = "Pool Size",
    y = "Number of Variants"
  )

pA

###############################################################################
# 5. FIGURE B — Ts/Tv RATIO BY POOL SIZE
###############################################################################

pB <- vars %>%
  ggplot(aes(x = Pool, y = TsTv, fill = Pool)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_classic(base_size = 14) +
  labs(
    title = "Transition/Transversion Ratio by Pool Size",
    x = "Pool Size",
    y = "Ts/Tv Ratio"
  )

pB
###############################################################################
# 6. FIGURE C — PROPORTION OF GENOME ≥ 10× BY POOL SIZE
###############################################################################

pC <- qc %>%
  ggplot(aes(x = Pool, y = Pct10x, fill = Pool)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic(base_size = 14) +
  labs(
    title = "Genome Coverage ≥10× by Pool Size",
    x = "Pool Size",
    y = "% Genome ≥10×"
  )
pC

###############################################################################
# 7. FIGURE D — MEAN DEPTH DISTRIBUTION BY POOL SIZE
###############################################################################

pD <- qc %>%
  ggplot(aes(x = Pool, y = MeanDepth, fill = Pool)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.5) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic(base_size = 14) +
  labs(
    title = "Mean Depth Across Pools",
    x = "Pool Size",
    y = "Mean Depth"
  )
pD

###############################################################################
# 8. COMPOSITE PANEL — FINAL PUBLICATION FIGURE
###############################################################################

final_fig <- ggarrange(
  pA, pB, pC, pD,
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2,
  common.legend = FALSE,
  font.label = list(size = 16)
)

final_fig

###############################################################################
# 9. EXPORT (optional)
###############################################################################

ggsave(
  "~/VSP-Training/Manuscript/BMJ-Sumission/Final_Figure_Pools_QC.png",
  final_fig,
  width = 12,
  height = 10,
  dpi = 300
)

###############################################################################
# END OF PIPELINE
###############################################################################
