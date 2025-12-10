###############################################################################
# FULL VISUAL QC PIPELINE FOR MANUSCRIPT FIGURES
# Author: Julien A. Nguinkal
# Goal: Produce clean, professional plots for QC and variant metrics,
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



















###############################################################################
# FULL VISUAL QC PIPELINE FOR MANUSCRIPT FIGURES
# Author: Julien A. Nguinkal
# Goal: Produce clean, professional plots for QC and variant metrics,
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
# 3b. COMMON THEME (NO LEGEND, CENTERED BOLD TITLES, BOLD X AXIS)
###############################################################################

theme_manuscript <- theme_classic(base_size = 14) +
  theme(
    legend.position = "none",                        # remove legend
    plot.title = element_text(
      hjust = 0, face = "bold"                    # centered, bold
    ),
    axis.title.x = element_text(face = "bold"),      # bold x-axis title
    axis.text.x  = element_text(face = "bold"), # bold x-axis tick labels
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")      # keep y normal (or set to bold if you prefer)
  )

###############################################################################
# 4. FIGURE A — DISTRIBUTION OF VARIANT COUNT BY POOL SIZE
###############################################################################

pA <- vars %>%
  ggplot(aes(x = Pool, y = TotalVariants, fill = Pool)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, color = "black", alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Total Variant Count",
    x = "Pool Size",
    y = "Number of Variants"
  ) +
  theme_manuscript

pA

###############################################################################
# 5. FIGURE B — Ts/Tv RATIO BY POOL SIZE
###############################################################################

pB <- vars %>%
  ggplot(aes(x = Pool, y = TsTv, fill = Pool)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    title = "Transition/Transversion",
    x = "Pool Size",
    y = "Ts/Tv Ratio"
  ) +
  theme_manuscript

pB

###############################################################################
# 6. FIGURE C — PROPORTION OF GENOME ≥ 10× BY POOL SIZE
###############################################################################

pC <- qc %>%
  ggplot(aes(x = Pool, y = Pct10x, fill = Pool)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "Genome Coverage ≥10×",
    x = "Pool Size",
    y = "% Genome ≥10×"
  ) +
  theme_manuscript

pC

###############################################################################
# 7. FIGURE D — MEAN DEPTH DISTRIBUTION BY POOL SIZE
###############################################################################

pD <- qc %>%
  ggplot(aes(x = Pool, y = MeanDepth, fill = Pool)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.5) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Mean Depth Across Pools",
    x = "Pool Size",
    y = "Mean Depth"
  ) +
  theme_manuscript

pD

###############################################################################
# 8. COMPOSITE PANEL — FINAL PUBLICATION FIGURE
###############################################################################

final_fig <- ggarrange(
  pA, pB, pC, pD,
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2,
  common.legend = FALSE,
  font.label = list(size = 16, face = "bold")
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

###############################################################################
# NATURE-STYLE QUALITY PANEL — NO BOXES / NO VIOLINS / STRIPLESS FACETS
###############################################################################

library(tidyverse)
library(janitor)

qc_qual <- qc_raw %>%
  clean_names() %>%
  filter(!is.na(pool_size)) %>%
  mutate(Pool = factor(pool_size, levels = sort(unique(pool_size))))

qc_qual_long <- qc_qual %>%
  transmute(
    Pool,
    `Q20 (%)`        = q20,
    `Q30 (%)`        = q30,
    `Mean quality`   = mean_qual,
    `Median quality` = median_qual
  ) %>%
  pivot_longer(
    cols      = -Pool,
    names_to  = "Metric",
    values_to = "Value"
  )

pQ_clean <- qc_qual_long %>%
  ggplot(aes(x = Pool, y = Value)) +
  
  # POINTS ONLY (jitter)
  geom_jitter(width = 0.10, size = 2.2, alpha = 0.7, color = "black") +
  
  # OPTIONAL: add mean line per facet
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    size = 0.4,
    color = "red",
    alpha = 0.7
  ) +
  
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  
  labs(
    title = "",
    x = "Pool Size",
    y = ""
  ) +
  
  theme_classic(base_size = 15) +
  
  # REMOVE facet strip background → very Nature-style
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 15, hjust = 0),
    
    # AXES
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(face = "bold", size = 13),
    axis.text.y  = element_text(size = 12),
    
    # REMOVE gridlines completely
    panel.grid = element_blank(),
    
    legend.position = "none",
    
    plot.title = element_text(hjust = 0, face = "bold")
  )

pQ_clean


# Export if needed
ggsave(
  "~/VSP-Training/Manuscript/BMJ-Sumission/Figure_Quality_Metrics_Stripless.png",
  pQ_clean,
  width = 10,
  height = 8,
  dpi = 300
)


# Unified colour palette for all figures
pool_colors <- c(
  "3"  = "#1b9e77",   # green
  "6"  = "#d95f02",   # orange
  "9"  = "#7570b3",   # purple
  "15" = "#e7298a"    # pink
)



###############################################################################
# NATURE-STYLE STRIPLESS PANEL — COLOURED BY POOL SIZE
###############################################################################

qc_qual <- qc_raw %>%
  clean_names() %>%
  filter(!is.na(pool_size)) %>%
  mutate(Pool = factor(pool_size, levels = sort(unique(pool_size))))

qc_qual_long <- qc_qual %>%
  transmute(
    Pool,
    `Q20 (%)`        = q20,
    `Q30 (%)`        = q30,
    `Mean quality`   = mean_qual,
    `Median quality` = median_qual
  ) %>%
  pivot_longer(
    cols      = -Pool,
    names_to  = "Metric",
    values_to = "Value"
  )

pQ_clean <- qc_qual_long %>%
  ggplot(aes(x = Pool, y = Value, color = Pool)) +
  
  # Points only (colored by pool)
  geom_jitter(width = 0.10, size = 2.8, alpha = 0.75) +
  
  # Optional mean bar (also colored by pool)
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    size = 0.45,
    aes(color = Pool),
    alpha = 0.8
  ) +
  
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  scale_color_manual(values = pool_colors) +
  
  labs(
    title = "",
    x = "Pool Size",
    y = ""
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 15, hjust = 0),
    
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(face = "bold", size = 13),
    axis.text.y  = element_text(size = 12),
    
    panel.grid = element_blank(),
    legend.position = "none",          # Pool color already encoded in points
    
    plot.title = element_text(hjust = 0, face = "bold")
  )

pQ_clean



###############################################################################
# QUALITY METRICS PANEL — NATURE-STYLE / STRIPLESS / NO BOXES
# Author: Dr. Julien A. Nguinkal
# Goal: Produce a clean, publication-ready panel comparing
#       Q20(%), Q30(%), mean quality, and median quality
#       across MPXV pooling levels (3, 6, 9, 15).
###############################################################################

library(tidyverse)
library(janitor)

set.seed(123)   # reproducible jitter

###############################################################################
# 1. LOAD RAW QC DATA
###############################################################################

qc_raw <- read.delim("~/VSP-Training/Manuscript/BMJ-Sumission/data_qc.csv")

###############################################################################
# 2. UNIFIED POOL COLOR PALETTE (used across all manuscript figures)
###############################################################################

pool_colors <- c(
  "3"  = "#1b9e77",   # green
  "6"  = "#d95f02",   # orange
  "9"  = "#7570b3",   # purple
  "15" = "#e7298a"    # pink
)

###############################################################################
# 3. CLEAN QC DATA
###############################################################################

qc_qual <- qc_raw %>%
  clean_names() %>%
  filter(!is.na(pool_size)) %>%
  mutate(
    Pool = factor(pool_size, levels = sort(unique(pool_size)))
  )

###############################################################################
# 4. CONVERT TO LONG FORMAT FOR FACETING
###############################################################################

qc_qual_long <- qc_qual %>%
  transmute(
    Pool,
    `Q20 (%)`        = q20,
    `Q30 (%)`        = q30,
    `Mean quality`   = mean_qual,
    `Median quality` = median_qual
  ) %>%
  pivot_longer(
    cols      = -Pool,
    names_to  = "Metric",
    values_to = "Value"
  )

###############################################################################
# 5. FINAL NATURE-STYLE PANEL (NO BOXES / FACETED / STRIPLESS)
###############################################################################

pQ_clean <- qc_qual_long %>%
  ggplot(aes(x = Pool, y = Value, color = Pool)) +
  
  geom_jitter(width = 0.10, size = 2.8, alpha = 0.75) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.45,
    size = 0.45,
    aes(color = Pool),
    alpha = 0.8
  ) +
  
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  
  scale_color_manual(values = pool_colors) +
  
  # FIXED HERE: removed na.translate
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0.05, 0.10))
  ) +
  
  labs(
    title = "",
    x = "Pool Size",
    y = ""
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 15, hjust = 0),
    
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(face = "bold", size = 13),
    axis.text.y  = element_text(size = 12),
    
    panel.grid = element_blank(),
    legend.position = "none"
  )


###############################################################################
# 6. DISPLAY FIGURE
###############################################################################

pQ_clean

###############################################################################
# 7. EXPORT
###############################################################################

ggsave(
  "~/VSP-Training/Manuscript/BMJ-Sumission/Figure_QualityMetrics_byPool.png",
  pQ_clean,
  width = 10,
  height = 8,
  dpi = 300
)

###############################################################################
# END OF SCRIPT
###############################################################################


###############################################################################
# FINAL NATURE-STYLE PANEL (NO BOXES, STRIPLESS, FREE Y AXES)
###############################################################################

pQ_clean <- qc_qual_long %>%
  ggplot(aes(x = Pool, y = Value, color = Pool)) +
  
  # Points only (no boxplot)
  geom_jitter(width = 0.10, size = 2.8, alpha = 0.75) +
  
  # Mean bar per pool
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.45,
    size = 0.45,
    aes(color = Pool),
    alpha = 0.8
  ) +
  
  # One facet per metric, each with its own y-scale
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  
  # Same pool colors as other figures
  scale_color_manual(values = pool_colors) +
  
  labs(
    title = "",
    x = "Pool Size",
    y = ""
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    # Stripless facet headers (Nature-style)
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 15, hjust = 0),
    
    # Axes
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(face = "bold", size = 13),
    axis.text.y  = element_text(size = 12),
    
    # Clean look
    panel.grid       = element_blank(),
    legend.position  = "none"
  )

pQ_clean







####################################

theme_manuscript <- theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      hjust = 0,
      face  = "bold"
    ),
    axis.title.x = element_text(face = "bold"),
    axis.text.x  = element_text(face = "bold"),
    axis.text.y  = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )



pA <- vars %>%
  ggplot(aes(x = Pool, y = TotalVariants, fill = Pool)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, color = "black", alpha = 0.5) +
  scale_fill_manual(values = pool_colors) +   # << here
  labs(
    title = "Total Variant Count",
    x = "Pool Size",
    y = "Number of Variants"
  ) +
  theme_manuscript
pA



pB <- vars %>%
  ggplot(aes(x = Pool, y = TsTv, fill = Pool)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6, color = "black") +
  scale_fill_manual(values = pool_colors) +   # << here
  labs(
    title = "Transition/Transversion",
    x = "Pool Size",
    y = "Ts/Tv Ratio"
  ) +
  theme_manuscript
pB


pC <- qc %>%
  ggplot(aes(x = Pool, y = Pct10x, fill = Pool)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.5, color = "black") +
  scale_fill_manual(values = pool_colors) +   # << here
  labs(
    title = "Genome Coverage ≥10×",
    x = "Pool Size",
    y = "% Genome ≥10×"
  ) +
  theme_manuscript
pC

pD <- qc %>%
  ggplot(aes(x = Pool, y = MeanDepth, fill = Pool)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.5, color = "black") +
  scale_fill_manual(values = pool_colors) +   # << here
  labs(
    title = "Mean Depth Across Pools",
    x = "Pool Size",
    y = "Mean Depth"
  ) +
  theme_manuscript
pD


final_fig <- ggarrange(
  pA, pB, pC, pD,
  ncol = 2, nrow = 2,
  common.legend = FALSE
)

final_fig





pQ_clean <- qc_qual_long %>%
  ggplot(aes(x = Pool, y = Value, color = Pool)) +
  
  # Points only
  geom_jitter(width = 0.10, size = 2.8, alpha = 0.75) +
  
  # Mean crossbar
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.45,
    size = 0.45,
    aes(color = Pool),
    alpha = 0.8
  ) +
  
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  
  scale_color_manual(values = pool_colors) +
  
  labs(
    title = "",
    x = "Pool Size",
    y = "Quality metric value"      # <── GLOBAL Y-AXIS TITLE ADDED
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 15, hjust = 0),
    
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),  # <── BOLD GLOBAL Y TITLE
    axis.text.x  = element_text(face = "bold", size = 13),
    axis.text.y  = element_text(size = 12),
    
    panel.grid      = element_blank(),
    legend.position = "none"
  )

pQ_clean






###############################################################################
# BUILD 4 SEPARATE PLOTS (EACH WITH ITS OWN Y-LABEL)
###############################################################################
qc_raw <- read.delim("~/VSP-Training/Manuscript/BMJ-Sumission/data_qc.csv")

qc <- qc_raw %>%
  clean_names() %>%
  filter(!is.na(pool_size)) %>%
  mutate(Pool = factor(pool_size, levels = sort(unique(pool_size))))
# Mean quality
p_mean <- qc %>%
  ggplot(aes(x = Pool, y = mean_qual, color = Pool)) +
  geom_jitter(width = 0.10, size = 3, alpha = 0.75) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.45, size = 0.5, aes(color = Pool)) +
  scale_color_manual(values = pool_colors) +
  labs(x = "Pool Size", y = "Mean quality") +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x  = element_text(face = "bold"),
    axis.text.y  = element_text(face = "bold")
  )

# Median quality
p_median <- qc %>%
  ggplot(aes(x = Pool, y = median_qual, color = Pool)) +
  geom_jitter(width = 0.10, size = 3, alpha = 0.75) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.45, size = 0.5, aes(color = Pool)) +
  scale_color_manual(values = pool_colors) +
  labs(x = "Pool Size", y = "Median quality") +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x  = element_text(face = "bold"),
    axis.text.y  = element_text(face = "bold")
  )

# Q20 (%)
p_q20 <- qc %>%
  ggplot(aes(x = Pool, y = q20, color = Pool)) +
  geom_jitter(width = 0.10, size = 3, alpha = 0.75) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.45, size = 0.5, aes(color = Pool)) +
  scale_color_manual(values = pool_colors) +
  labs(x = "Pool Size", y = "Q20 (%)") +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x  = element_text(face = "bold"),
    axis.text.y  = element_text(face = "bold")
  )

# Q30 (%)
p_q30 <- qc %>%
  ggplot(aes(x = Pool, y = q30, color = Pool)) +
  geom_jitter(width = 0.10, size = 3, alpha = 0.75) +
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.45, size = 0.5, aes(color = Pool)) +
  scale_color_manual(values = pool_colors) +
  labs(x = "Pool Size", y = "Q30 (%)") +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x  = element_text(face = "bold"),
    axis.text.y  = element_text(face = "bold")
  )

###############################################################################
# COMBINE IN 2×2 PANEL (Nature-style)
###############################################################################

final_quality_panel <- (p_mean | p_median) /
  (p_q20  | p_q30)

final_quality_panel



library(tidyverse)
library(janitor)
library(ggrepel)

###############################################################################
# 1. Load and clean QC data
###############################################################################

qc_raw <- read.delim("~/VSP-Training/Manuscript/BMJ-Sumission/data_qc.csv")

qc <- qc_raw %>%
  clean_names() %>% 
  # expected columns: sample, pool_size, mean_depth, genome_ge_10x_pct,
  #                   q20, q30, mean_qual, median_qual, ...
  filter(!is.na(pool_size), !is.na(sample)) %>%
  mutate(
    PoolSize = as.numeric(pool_size)
  )

###############################################################################
# 2. Restrict to samples present in ≥ 2 pools (the informative ones)
###############################################################################

multi_pool_samples <- qc %>%
  group_by(sample) %>%
  summarise(n_pools = n_distinct(PoolSize), .groups = "drop") %>%
  filter(n_pools >= 2) %>%
  pull(sample)

qc_multi <- qc %>%
  filter(sample %in% multi_pool_samples)

# Inspect which samples are used:
unique(qc_multi$sample)
# e.g. "MPXV-S03", "MPXV-S04", "MPXV-S05", "MPXV-S06", "MPXV-S07", "MPXV-S19", ...

###############################################################################
# 3. Long format for the 4 metrics of interest
###############################################################################

qc_long <- qc_multi %>%
  select(
    sample, PoolSize,
    mean_depth,
    genome_ge_10x_pct,
    q30,
    median_qual
  ) %>%
  pivot_longer(
    cols = c(mean_depth, genome_ge_10x_pct, q30, median_qual),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("mean_depth", "genome_ge_10x_pct", "q30", "median_qual"),
      labels = c(
        "Mean depth",
        "Genome ≥10× (%)",
        "Q30 (%)",
        "Median quality"
      )
    )
  )

###############################################################################
# 4. Base plot: trend per sample across pool sizes, facetted by metric
###############################################################################

# Optionally highlight one or two key samples more strongly
highlight_samples <- c("MPXV-S03", "MPXV-S05")  # adjust as you like

p_trend <- ggplot(qc_long, aes(x = PoolSize, y = Value, group = sample)) +
  
  # All samples in light grey (context)
  geom_line(color = "grey80", size = 0.7, alpha = 0.8) +
  geom_point(color = "grey70", size = 2, alpha = 0.8) +
  
  # Highlighted samples in colour
  geom_line(
    data = subset(qc_long, sample %in% highlight_samples),
    aes(color = sample),
    size = 1.1,
    alpha = 0.95
  ) +
  geom_point(
    data = subset(qc_long, sample %in% highlight_samples),
    aes(color = sample),
    size = 2.5,
    alpha = 0.95
  ) +
  
  # Optional labels at the largest pool size for highlighted samples
  geom_label_repel(
    data = qc_long %>%
      filter(sample %in% highlight_samples) %>%
      group_by(sample, Metric) %>%
      filter(PoolSize == max(PoolSize)) %>%
      ungroup(),
    aes(label = sample, color = sample),
    size = 3,
    label.size = 0.1,
    label.r = unit(0.15, "lines"),
    segment.size = 0.3,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  
  scale_x_continuous(
    breaks = c(3, 6, 9, 15),
    minor_breaks = NULL
  ) +
  
  labs(
    x = "Pool size",
    y = "Value per sample"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 14, hjust = 0),
    axis.title.x     = element_text(face = "bold"),
    axis.title.y     = element_text(face = "bold"),
    axis.text.x      = element_text(face = "bold"),
    legend.position  = "none"
  )

p_trend






library(tidyverse)
library(janitor)
library(ggrepel)

###############################################################################
# 1. Load and clean QC data
###############################################################################

qc_raw <- read.delim("~/VSP-Training/Manuscript/BMJ-Sumission/data_qc.csv")

qc <- qc_raw %>%
  clean_names() %>%
  filter(!is.na(pool_size), !is.na(sample)) %>%
  mutate(PoolSize = as.numeric(pool_size))

###############################################################################
# 2. Keep only samples present in 2 or 3 pools
###############################################################################

multi_pool_samples <- qc %>%
  group_by(sample) %>%
  summarise(n_pools = n_distinct(PoolSize), .groups = "drop") %>%
  filter(n_pools >= 2, n_pools <= 3) %>%    # exactly 2 or 3 pools
  pull(sample)

qc_multi <- qc %>%
  filter(sample %in% multi_pool_samples)

###############################################################################
# 3. Long format for the 4 metrics of interest
###############################################################################

qc_long <- qc_multi %>%
  select(
    sample, PoolSize,
    mean_depth,
    genome_ge_10x_pct,
    q30,
    median_qual
  ) %>%
  pivot_longer(
    cols = c(mean_depth, genome_ge_10x_pct, q30, median_qual),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("mean_depth", "genome_ge_10x_pct", "q30", "median_qual"),
      labels = c(
        "Mean depth",
        "Genome ≥10× (%)",
        "Q30 (%)",
        "Median quality"
      )
    )
  )

###############################################################################
# 4. Figure 3 – Only shared samples (2–3 pools), lines coloured by sample
###############################################################################

p_fig3 <- ggplot(qc_long, aes(x = PoolSize, y = Value, group = sample, color = sample)) +
  geom_line(size = 1.0, alpha = 0.9) +
  geom_point(size = 2.4, alpha = 0.9) +
  
  # label each sample at its largest pool size
  geom_label_repel(
    data = qc_long %>%
      group_by(sample, Metric) %>%
      filter(PoolSize == max(PoolSize)) %>%
      ungroup(),
    aes(label = sample),
    size = 3,
    label.size = 0.15,
    label.r = unit(0.15, "lines"),
    segment.size = 0.3,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  
  scale_x_continuous(
    breaks = c(3, 6, 9, 15),
    minor_breaks = NULL
  ) +
  
  labs(
    x = "Pool size",
    y = "Value per sample"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 14, hjust = 0),
    axis.title.x     = element_text(face = "bold"),
    axis.title.y     = element_text(face = "bold"),
    axis.text.x      = element_text(face = "bold"),
    legend.position  = "none"
  )

p_fig3








library(tidyverse)
library(janitor)
library(ggrepel)
library(patchwork)

###############################################################################
# 1. Load and clean QC data
###############################################################################

qc_raw <- read.delim("~/VSP-Training/Manuscript/BMJ-Sumission/data_qc.csv")

qc <- qc_raw %>%
  clean_names() %>%
  filter(!is.na(pool_size), !is.na(sample)) %>%
  mutate(
    PoolSize = as.numeric(pool_size)
  )

###############################################################################
# 2. Keep only samples present in ≥ 2 pools
###############################################################################

multi_pool_samples <- qc %>%
  group_by(sample) %>%
  summarise(n_pools = n_distinct(PoolSize), .groups = "drop") %>%
  filter(n_pools >= 2) %>%
  pull(sample)

qc_multi <- qc %>%
  filter(sample %in% multi_pool_samples)

###############################################################################
# 3. Prepare long format for metrics of interest
###############################################################################

qc_long <- qc_multi %>%
  select(
    sample, PoolSize,
    mean_depth,
    genome_ge_10x_pct,
    q30,
    median_qual
  ) %>%
  pivot_longer(
    cols = c(mean_depth, genome_ge_10x_pct, q30, median_qual),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("mean_depth", "genome_ge_10x_pct", "q30", "median_qual"),
      labels = c(
        "Mean depth",
        "Genome ≥10× (%)",
        "Q30 (%)",
        "Median quality"
      )
    )
  )

###############################################################################
# 4. Helper: make one plot per metric (with its own y-axis label)
###############################################################################

make_trend_plot <- function(df_metric, y_lab) {
  ggplot(df_metric, aes(x = PoolSize, y = Value, group = sample, color = sample)) +
    geom_line(size = 1.0, alpha = 0.9) +
    geom_point(size = 2.4, alpha = 0.9) +
    geom_label_repel(
      data = df_metric %>%
        group_by(sample) %>%
        filter(PoolSize == max(PoolSize)) %>%
        ungroup(),
      aes(label = sample),
      size         = 3,
      label.size   = 0.15,
      label.r      = unit(0.15, "lines"),
      segment.size = 0.3,
      show.legend  = FALSE,
      max.overlaps = Inf
    ) +
    scale_x_continuous(
      breaks = c(3, 6, 9, 15),
      minor_breaks = NULL
    ) +
    labs(
      x = "Pool size",
      y = y_lab
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title.x    = element_text(face = "bold"),
      axis.title.y    = element_text(face = "bold"),
      axis.text.x     = element_text(face = "bold"),
      axis.text.y     = element_text(face = "bold")
    )
}

###############################################################################
# 5. Build the 4 subplots
###############################################################################

p_mean_depth <- qc_long %>%
  filter(Metric == "Mean depth") %>%
  make_trend_plot(y_lab = "Mean depth")

p_cov10x <- qc_long %>%
  filter(Metric == "Genome ≥10× (%)") %>%
  make_trend_plot(y_lab = "Genome ≥10× (%)")

p_q30 <- qc_long %>%
  filter(Metric == "Q30 (%)") %>%
  make_trend_plot(y_lab = "Q30 (%)")

p_med_qual <- qc_long %>%
  filter(Metric == "Median quality") %>%
  make_trend_plot(y_lab = "Median quality")

###############################################################################
# 6. Combine into Figure 3 (2×2 panel)
###############################################################################

figure3_trends <- (p_mean_depth | p_cov10x) /
  (p_q30        | p_med_qual)

figure3_trends

















###############################################################################
# Figure 3 – Three alternative visualizations of pooling impact
# Data: data_qc.csv  (sample, pool_size, mean_depth, genome_ge_10x_pct, q30, mean_qual, median_qual)
###############################################################################

library(tidyverse)
library(janitor)
library(ggrepel)
library(patchwork)

###############################################################################
# 1. LOAD & CLEAN DATA
###############################################################################

qc_raw <- read.delim("~/VSP-Training/Manuscript/BMJ-Sumission/data_qc.csv")

qc <- qc_raw %>%
  clean_names() %>%                              # sample, pool_size, mean_depth, genome_ge_10x_pct, q30, mean_qual, median_qual
  filter(!is.na(pool_size), !is.na(sample)) %>%
  mutate(
    PoolSize = as.numeric(pool_size)
  )

# Samples present in ≥ 2 pools (the ones of interest for dilution trends)
multi_pool_samples <- qc %>%
  group_by(sample) %>%
  summarise(n_pools = n_distinct(PoolSize), .groups = "drop") %>%
  filter(n_pools >= 2) %>%
  pull(sample)

qc_multi <- qc %>%
  filter(sample %in% multi_pool_samples)

# Long format for the 4 key metrics
qc_long <- qc_multi %>%
  select(
    sample, PoolSize,
    mean_depth,
    genome_ge_10x_pct,
    q30,
    median_qual
  ) %>%
  pivot_longer(
    cols = c(mean_depth, genome_ge_10x_pct, q30, median_qual),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("mean_depth", "genome_ge_10x_pct", "q30", "median_qual"),
      labels = c(
        "Mean depth",
        "Genome ≥10× (%)",
        "Q30 (%)",
        "Median quality"
      )
    )
  )

# A small color palette for samples (only repeated ones, so not huge)
sample_colors <- scales::hue_pal()(length(unique(qc_long$sample)))
names(sample_colors) <- sort(unique(qc_long$sample))

###############################################################################
# OPTION 1 – PURE SLOPE PLOTS (ONE LINE PER SAMPLE, FACETTED BY METRIC)
###############################################################################

p_opt1 <- ggplot(qc_long, aes(x = PoolSize, y = Value, group = sample, color = sample)) +
  geom_line(size = 1.0, alpha = 0.9) +
  geom_point(size = 2.4, alpha = 0.9) +
  geom_label_repel(
    data = qc_long %>%
      group_by(sample, Metric) %>%
      filter(PoolSize == max(PoolSize)) %>%
      ungroup(),
    aes(label = sample),
    size         = 3,
    label.size   = 0.15,
    label.r      = unit(0.15, "lines"),
    segment.size = 0.3,
    show.legend  = FALSE,
    max.overlaps = Inf
  ) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  scale_x_continuous(
    breaks = c(3, 6, 9, 15),
    minor_breaks = NULL
  ) +
  scale_color_manual(values = sample_colors) +
  labs(
    x = "Pool size",
    y = "Value per sample"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 14, hjust = 0),
    axis.title.x     = element_text(face = "bold"),
    axis.title.y     = element_text(face = "bold"),
    axis.text.x      = element_text(face = "bold"),
    axis.text.y      = element_text(face = "bold")
  )

p_opt1    # → Option 1: direct trajectories per sample across pool sizes


###############################################################################
# OPTION 2 – RELATIVE CHANGE VS SMALLEST POOL (NORMALISED PER SAMPLE & METRIC)
###############################################################################

qc_rel <- qc_long %>%
  group_by(sample, Metric) %>%
  mutate(
    base_pool = PoolSize[which.min(PoolSize)],
    base_val  = Value[which.min(PoolSize)],
    RelValue  = ifelse(base_val > 0, Value / base_val, NA_real_)
  ) %>%
  ungroup()

p_opt2 <- ggplot(qc_rel, aes(x = PoolSize, y = RelValue, group = sample, color = sample)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", size = 0.4) +
  geom_line(size = 1.0, alpha = 0.9) +
  geom_point(size = 2.4, alpha = 0.9) +
  geom_label_repel(
    data = qc_rel %>%
      group_by(sample, Metric) %>%
      filter(PoolSize == max(PoolSize)) %>%
      ungroup(),
    aes(label = sample),
    size         = 3,
    label.size   = 0.15,
    label.r      = unit(0.15, "lines"),
    segment.size = 0.3,
    show.legend  = FALSE,
    max.overlaps = Inf
  ) +
  facet_wrap(~ Metric, ncol = 2) +
  scale_x_continuous(
    breaks = c(3, 6, 9, 15),
    minor_breaks = NULL
  ) +
  scale_y_continuous(
    breaks = seq(0, 1.6, 0.2),
    limits = c(0, NA)
  ) +
  scale_color_manual(values = sample_colors) +
  labs(
    x = "Pool size",
    y = "Relative value (vs smallest pool)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 14, hjust = 0),
    axis.title.x     = element_text(face = "bold"),
    axis.title.y     = element_text(face = "bold"),
    axis.text.x      = element_text(face = "bold"),
    axis.text.y      = element_text(face = "bold")
  )

p_opt2    # → Option 2: relative curves (1.0 = no change vs smallest pool)


###############################################################################
# OPTION 3 – BOX PLOTS BY POOL + HIGHLIGHTED REPEATED SAMPLES
###############################################################################

# For this, we use ALL samples (qc), and overlay only the repeated ones (qc_multi).
qc_long_all <- qc %>%
  select(
    sample, PoolSize,
    mean_depth,
    genome_ge_10x_pct,
    q30,
    median_qual
  ) %>%
  pivot_longer(
    cols = c(mean_depth, genome_ge_10x_pct, q30, median_qual),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("mean_depth", "genome_ge_10x_pct", "q30", "median_qual"),
      labels = c(
        "Mean depth",
        "Genome ≥10× (%)",
        "Q30 (%)",
        "Median quality"
      )
    )
  )

qc_long_multi <- qc_long %>%
  mutate(sample = factor(sample))

p_opt3 <- ggplot() +
  # Background distribution (all samples) as boxplots per pool
  geom_boxplot(
    data = qc_long_all,
    aes(x = factor(PoolSize), y = Value),
    fill = "grey90",
    color = "black",
    width = 0.55,
    outlier.shape = NA
  ) +
  # Highlight repeated samples as coloured points with connecting lines
  geom_line(
    data = qc_long_multi,
    aes(x = PoolSize, y = Value, group = sample, color = sample),
    size = 0.9,
    alpha = 0.9
  ) +
  geom_point(
    data = qc_long_multi,
    aes(x = PoolSize, y = Value, color = sample),
    size = 2.6,
    alpha = 0.95
  ) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  scale_color_manual(values = sample_colors) +
  labs(
    x = "Pool size",
    y = "Value"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 14, hjust = 0),
    axis.title.x     = element_text(face = "bold"),
    axis.title.y     = element_text(face = "bold"),
    axis.text.x      = element_text(face = "bold"),
    axis.text.y      = element_text(face = "bold")
  )

p_opt3    # → Option 3: intuitive boxplots + highlighted repeated samples


###############################################################################
# OPTIONAL – ARRANGE SIDE BY SIDE TO COMPARE THE THREE OPTIONS
###############################################################################

# For example:
(p_opt1 | p_opt2) / p_opt3







###############################################################################
# FIGURE 3 – PAIRED SLOPE PLOTS (ONLY 2 SAMPLES PER METRIC)
# - Uses data_qc.csv
# - For each metric, keeps at most:
#     1 sample present in ≥ 3 pools
#     1 sample present in exactly 2 pools
# - One plot per metric, each with its own x / y axis labels
###############################################################################

library(tidyverse)
library(janitor)
library(ggrepel)

###############################################################################
# 1. Load & clean data
###############################################################################

qc_raw <- read.delim("~/VSP-Training/Manuscript/BMJ-Sumission/data_qc.csv")

qc <- qc_raw %>%
  clean_names() %>%    # sample, pool_size, mean_depth, genome_ge_10x_pct, q30, mean_qual, median_qual
  filter(!is.na(pool_size), !is.na(sample)) %>%
  mutate(
    PoolSize = as.numeric(pool_size),
    sample   = as.character(sample)
  )

###############################################################################
# 2. Keep only samples present in ≥ 2 pools (for dilution trends)
###############################################################################

multi_pool_samples <- qc %>%
  group_by(sample) %>%
  summarise(n_pools = n_distinct(PoolSize), .groups = "drop") %>%
  filter(n_pools >= 2) %>%
  pull(sample)

qc_multi <- qc %>%
  filter(sample %in% multi_pool_samples)

###############################################################################
# 3. Long format for the 4 key metrics
###############################################################################

qc_long <- qc_multi %>%
  select(
    sample, PoolSize,
    mean_depth,
    genome_ge_10x_pct,
    q30,
    median_qual
  ) %>%
  pivot_longer(
    cols = c(mean_depth, genome_ge_10x_pct, q30, median_qual),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c("mean_depth", "genome_ge_10x_pct", "q30", "median_qual"),
      labels = c(
        "Mean depth",
        "Genome ≥10× (%)",
        "Q30 (%)",
        "Median quality"
      )
    )
  )

# Color palette for samples (consistent across all plots)
sample_colors <- scales::hue_pal()(length(unique(qc_long$sample)))
names(sample_colors) <- sort(unique(qc_long$sample))

###############################################################################
# 4. Helper: pick at most 2 samples per metric (1 with ≥3 pools, 1 with 2 pools)
###############################################################################

select_two_samples_for_metric <- function(df_metric) {
  # df_metric: subset of qc_long for a single Metric
  summary_tbl <- df_metric %>%
    group_by(sample) %>%
    summarise(
      n_pools = n_distinct(PoolSize),
      mean_val = mean(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # First: sample present in ≥ 3 pools (if any), choose the one with highest mean_val
  cand3 <- summary_tbl %>%
    filter(n_pools >= 3) %>%
    arrange(desc(mean_val)) %>%
    slice_head(n = 1)
  
  # Second: sample present in exactly 2 pools, not already chosen
  cand2 <- summary_tbl %>%
    filter(n_pools == 2, !(sample %in% cand3$sample)) %>%
    arrange(desc(mean_val)) %>%
    slice_head(n = 1)
  
  chosen_samples <- bind_rows(cand3, cand2) %>%
    pull(sample) %>%
    unique()
  
  df_metric %>%
    filter(sample %in% chosen_samples)
}

###############################################################################
# 5. Helper: make one clean slope plot for a given metric
###############################################################################

make_metric_plot <- function(df_metric, y_lab) {
  ggplot(df_metric, aes(x = PoolSize, y = Value, group = sample, color = sample)) +
    geom_line(size = 1.1, alpha = 0.9) +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_label_repel(
      data = df_metric %>%
        group_by(sample) %>%
        filter(PoolSize == max(PoolSize)) %>%
        ungroup(),
      aes(label = sample),
      size         = 3,
      label.size   = 0.2,
      label.r      = unit(0.15, "lines"),
      segment.size = 0.3,
      show.legend  = FALSE,
      max.overlaps = Inf
    ) +
    scale_x_continuous(
      breaks = c(3, 6, 9, 15),
      minor_breaks = NULL
    ) +
    scale_color_manual(values = sample_colors) +
    labs(
      x = "Pool size",
      y = y_lab
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title.x    = element_text(face = "bold"),
      axis.title.y    = element_text(face = "bold"),
      axis.text.x     = element_text(face = "bold"),
      axis.text.y     = element_text(face = "bold")
    )
}

###############################################################################
# 6. Build the 4 subplots (each with only 2 samples)
###############################################################################

## Mean depth
df_mean_depth <- qc_long %>%
  filter(Metric == "Mean depth") %>%
  select_two_samples_for_metric()

p_mean_depth <- make_metric_plot(df_mean_depth, y_lab = "Mean depth")

## Genome ≥10× (%)
df_cov10x <- qc_long %>%
  filter(Metric == "Genome ≥10× (%)") %>%
  select_two_samples_for_metric()

p_cov10x <- make_metric_plot(df_cov10x, y_lab = "Genome ≥10× (%)")

## Q30 (%)
df_q30 <- qc_long %>%
  filter(Metric == "Q30 (%)") %>%
  select_two_samples_for_metric()

p_q30 <- make_metric_plot(df_q30, y_lab = "Q30 (%)")

## Median quality
df_med <- qc_long %>%
  filter(Metric == "Median quality") %>%
  select_two_samples_for_metric()

p_med_qual <- make_metric_plot(df_med, y_lab = "Median quality")

###############################################################################
# 7. Print the 4 panels (you can arrange them later as A/B/C/D)
###############################################################################

p_mean_depth
p_cov10x
p_q30
p_med_qual

library(patchwork)

panel_F3 <- (p_mean_depth | p_cov10x) /
  (p_q30        | p_med_qual)

panel_F3


