# Load libraries
library(ggplot2)
library(dplyr)

# Create data frame
qc_df <- data.frame(
  Run = factor(c("RUN-1", "RUN-2", "RUN-3", "RUN-4"), levels = c("RUN-1", "RUN-2", "RUN-3", "RUN-4")),
  SampleSize = c("n=3", "n=6", "n=9", "n=15"),
  Q20 = c(96.58, 98.18, 98.18, 97.79),
  Q30 = c(91.44, 94.52, 94.52, 94.53)
)

# Pivot to long format
qc_long <- qc_df %>%
  tidyr::pivot_longer(cols = c(Q20, Q30), names_to = "Quality", values_to = "Percentage")

# Combine run and sample size
qc_long$RunLabel <- paste0(qc_long$Run, "\n", qc_long$SampleSize)

# Define colors (Lancet/Nature-like)
colors <- c("Q20" = "#1f77b4", "Q30" = "#ff7f0e")

# Plot
ggplot(qc_long, aes(x = RunLabel, y = Percentage, fill = Quality)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = sprintf("%.2f%%", Percentage)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = colors) +
  labs(
    title = "Mean Q20 and Q30 Read Proportions per Sequencing Run",
    x = NULL,
    y = "Percentage of Reads",
    fill = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 14),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray70")
    
    
  )












# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Create data frame
qc_df <- data.frame(
  Run = factor(c("RUN-1", "RUN-2", "RUN-3", "RUN-4"), levels = c("RUN-1", "RUN-2", "RUN-3", "RUN-4")),
  SampleSize = c("n=3", "n=6", "n=9", "n=15"),
  Q20 = c(96.58, 98.18, 98.18, 97.79),
  Q30 = c(91.44, 95.52, 94.52, 94.53)
)

# Pivot to long format
qc_long <- qc_df %>%
  pivot_longer(cols = c(Q20, Q30), names_to = "Quality", values_to = "Percentage")

# Combine run and sample size
qc_long$RunLabel <- paste0(qc_long$Run, "\n", qc_long$SampleSize)

# Define custom colors
colors <- c("Q20" = "#1f77b4", "Q30" = "#ff7f0e")

# Create plot object
p <- ggplot(qc_long, aes(x = RunLabel, y = Percentage, fill = Quality)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = sprintf("%.2f%%", Percentage)),
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = colors) +
  labs(
    title = "Mean Q20 and Q30 Read Proportions per Sequencing Run",
    x = NULL,
    y = "Percentage of Reads",
    fill = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 14),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray70")
  )

# Print plot to screen
print(p)

# Optional: Save to high-quality PDF
ggsave("QScore_Distribution_by_Run.pdf", plot = p, width = 10, height = 6, dpi = 300)






# Load required libraries
library(ggplot2)
library(readr)
library(dplyr)

# Load the data (adjust path if necessary)
df <- read_csv("repeated_sample_assembly_metrics.csv")

# Ensure proper ordering of runs
df$Run <- factor(df$Run, levels = c("Run-03", "Run-06", "Run-09"))

# Optional: order metric levels explicitly
df$Metric <- factor(df$Metric, levels = c("N50", "Largest Contig"))

# Main plot
ggplot(df, aes(x = Run, y = Value_Kbp, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  facet_wrap(~ Sample, ncol = 3, labeller = labeller(Sample = function(x) paste("Sample-ID:", x))) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Run",
    y = "Assembly Metric (Kbp)",
    fill = "Metric"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(size = 0.3),
    panel.grid.minor = element_blank()
  )
