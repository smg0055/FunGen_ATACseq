# Load necessary packages
library(ggplot2)
library(stringr)
library(gridExtra)

# Import data
heatstress_enrichment_data <- read.csv("./heatstressknownResults.csv", header = TRUE, fill=TRUE)
control_enrichment_data <- read.csv("./controlknownResults.csv", header = TRUE, fill=TRUE)

# Select top 10 motifs
top_heatstress_enrichment_data <- head(heatstress_enrichment_data, 10)
top_control_enrichment_data <- head(control_enrichment_data, 10)

#shorten the names for labels
heatstress_short_names <- str_extract(top_heatstress_enrichment_data$Motif.Name, "[^/]+")
control_short_names <- str_extract(top_control_enrichment_data$Motif.Name, "[^/]+")

# Create bar plot
heatstress_plot <- ggplot(top_heatstress_enrichment_data, aes(y = reorder(heatstress_short_names, `Log.P.value`), x = `Log.P.value`, fill = `heatstress_short_names`)) +
  geom_bar(stat = "identity") +
  scale_x_reverse() +
  ylab("Motif") +
  xlab("Log P-value") +
  ggtitle("Top 10 Enriched Motifs by 
  Log P-value Under Heat Stress") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# Create bar plot
control_plot <- ggplot(top_control_enrichment_data, aes(y = reorder(control_short_names, `Log.P.value`), x = `Log.P.value`, fill = `control_short_names`)) +
  geom_bar(stat = "identity") +
  scale_x_reverse() +
  ylab("Motif") +
  xlab("Log P-value") +
  ggtitle("Top 10 Enriched Motifs by 
  Log P-value From Control") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# Arrange the plots side by side
arrange_plots <- grid.arrange(heatstress_plot, control_plot, ncol=2)

# Export the arranged plot
ggsave("motifenrichmentbargraphs.jpg", arrange_plots, width=24, height=10, units="cm")
