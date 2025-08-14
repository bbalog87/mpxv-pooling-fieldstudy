# Install necessary packages if not already installed
if (!require(ggplot2)) install.packages("ggplot2")

# Data for taxa distribution
taxa <- c("Monkeypox virus", "Orthopoxvirus", "Orthopoxvirus", "Orthopoxvirus", "Camelpox virus",
          "Orthopoxvirus", "Orthopoxvirus", "Orthopoxvirus", "Orthopoxvirus", "Orthopoxvirus",
          "Monkeypox virus", "Orthopoxvirus", "Orthopoxvirus", "Orthopoxvirus", "Monkeypox virus",
          "Monkeypox virus", "Orthopoxvirus", "Monkeypox virus", "Monkeypox virus", "Orthopoxvirus",
          "Monkeypox virus", "Orthopoxvirus", "Vaccinia virus", "Orthopoxvirus", "Orthopoxvirus", 
          "Orthopoxvirus", "Orthopoxvirus", "Monkeypox virus", "Orthopoxvirus", "Orthopoxvirus",
          "Orthopoxvirus", "Monkeypox virus", "Monkeypox virus", "Monkeypox virus", "Monkeypox virus",
          "Monkeypox virus", "Orthopoxvirus", "Monkeypox virus", "Monkeypox virus", "Monkeypox virus",
          "Monkeypox virus")

# Install the necessary packages if not already installed
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(ggsci)) install.packages("ggsci")
library(ggplot2)
library(dplyr)
library(ggsci)

# Create a data frame and filter rows with Freq > 1000
df <- as.data.frame(table(kaiju_taxa$V1)) %>% filter(Freq > 1000)

# Calculate the proportion of each taxon
df$Proportion <- df$Freq / sum(df$Freq) * 100  # Convert proportion to percentage

# Sort taxa by proportion in decreasing order
df <- df %>% arrange(desc(Proportion))

# Group all taxa after the 10th as "Other"
df$Var1 <- as.character(df$Var1)
df$Var1[11:nrow(df)] <- "Other"

# Recalculate proportions after grouping "Other"
df <- df %>%
  group_by(Var1) %>%
  summarise(Proportion = sum(Proportion)) %>%
  ungroup()

# Reorder taxa (Var1) by proportion in decreasing order to ensure both bar and legend match
df$Var1 <- factor(df$Var1, levels = df$Var1[order(df$Proportion, decreasing = FALSE)])

# Create a color palette: npg colors for top 10, dark grey for "Other"
palette_colors <- rev(c(ggsci::pal_npg("nrc")(10), "darkgrey"))

# Create a stacked bar chart with the y-axis as a percentage and sorted legend
ggplot(df, aes(x = "Taxa Distribution", y = Proportion, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(title = "Viral Taxa Distribution from Kaiju Results", x = "", y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 18),     # Bold y-axis title
    axis.text.y = element_text(size = 18),                    # Set y-axis text size
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  # Center and bold title
    legend.text = element_text(face = "bold", size = 18),      # Bold and larger legend text
    legend.title = element_blank()                             # Remove legend title
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +  # Set y-axis limits and breaks
  scale_fill_manual(values = palette_colors)



# Read the data
data <- allSamples_kaiju_results

# Name the columns for easier reference
colnames(data) <- c("SampleGroup", "Category", "Identifier", "TaxID", "Species")

# Calculate the proportion of each species within each sample group
summary <- data %>% 
  group_by(SampleGroup, Species) %>%
  summarise(Count = n()) %>% filter(Count >500) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  arrange(SampleGroup, desc(Proportion))

# Print the summary
print(summary)





# Load necessary libraries
library(dplyr)
library(stringr)


data <-allSamples_kaiju_contigs_results
# Rename columns for easier reference
colnames(data) <- c("SampleGroup", "Status", "ContigInfo", "TaxID", "Species")

# Extract contig lengths from the ContigInfo column
data <- data %>%
  mutate(ContigLength = as.numeric(str_extract(ContigInfo, "(?<=_length_)[0-9]+"))) %>%
  filter(ContigLength >= 200)  # Filter out contigs with length less than 200

# Calculate total contig length and species proportion for each sample group
summary <- data %>%
  group_by(SampleGroup) %>%
  summarise(TotalLength = sum(ContigLength)) %>%
  inner_join(
    data %>%
      group_by(SampleGroup, Species) %>%
      summarise(SpeciesLength = sum(ContigLength)) %>%
      mutate(Proportion = SpeciesLength / sum(SpeciesLength)),
    by = "SampleGroup"
  ) %>% arrange(SampleGroup, Proportion)

# Print the summary
print(summary)


coverage_all_runs<-combined_coverage
colnames(coverage_all_runs)<-c("Run", "Sample", "Ref", "Position", "Depth")



# Read the data file
#coverage_data <- read_tsv("coverage_all_runs.txt", col_names = c("Run", "Sample", "Reference", "Position", "Coverage"))

# Filter data for specific samples (e.g., S1, S7, S13)
coverage_data_filtered <- coverage_all_runs %>%
  filter(Sample %in% c("S7_mpxv", "S13_mpxv", "S114_mpxv")) %>%
  mutate(Run = factor(Run), Sample = factor(Sample))

# Plot depth across positions for each sample and run
ggplot(coverage_data_filtered,
       aes(x = Position, y = Depth, color = Run, group = Run)) +
  geom_line(alpha = 0.8, size = 1) +
  facet_wrap(~ Sample, scales = "free_x") +
  theme_minimal() +
  labs(title = "Comparative Depth Coverage for Samples S1, S7, and S13 Across Runs",
       x = "Genomic Position",
       y = "Coverage",
       color = "Run") +
  scale_color_manual(values = c("RUN1" = "#1f78b4", "RUN2" = "#33a02c", "RUN3" = "#e31a1c"))  # Custom colors for clear distinction

