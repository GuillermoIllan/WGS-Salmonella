#!/usr/bin/env Rscript

# ----------------------------
# Setup
# ----------------------------

setwd("/run/user/1000/gvfs/smb-share:server=bk-mxq.local,share=genomica/Proyectos_Guille/WGS_Salmonella")

# Load necessary libraries
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(plotly)
library(ggnewscale)
library(ggpubr)
library(grid)  # for unit()
# Note: ggsave uses ragg devices below, so ragg must be installed
# install.packages("ragg") if needed

# ----------------------------
# Load metadata
# ----------------------------

# Load metadata
metadata <- read.table("sistr/sistr-output.tab", header = TRUE, sep = "\t")

# Remove NS2 because it is not Salmonella
metadata <- metadata %>%
  filter(genome != "NS2")

# Harmonize label for the sample point
metadata$sample_point <- gsub(
  "Irrigation community water",
  "Reservoir water",
  metadata$sample_point
)

# Read the Newick tree file
tree <- read.tree("kSNP4/RAxML_bestTree.core_raxml.tree")

# ----------------------------
# Add metadata to the tree (fan tree + metadata rings)
# ----------------------------

# Add metadata to the tree
# tree$tip.label <- metadata$genome[match(tree$tip.label, metadata$genome)]
tree$tip.label

meta_filt <- metadata %>% filter(genome %in% tree$tip.label)
length(tree$tip.label) == nrow(meta_filt)

good_order <- match(tree$tip.label, meta_filt$genome)
meta_filt <- meta_filt[good_order, ]

meta_filt$serogroup <- as.factor(meta_filt$serogroup)
rownames(meta_filt) <- meta_filt$genome

p <- ggtree(tree, layout = "fan", branch.length = "none") +
  geom_tiplab(size = 5) +
  ggtitle("Core Genome SNPs Maximum-Likelihood Phylogenetic Tree")

# Ring 1: Subspecies
cc <- data.frame(as.character(meta_filt$cgmlst_subspecies), row.names = tree$tip.label)
colnames(cc) <- "Subspecies"
cc_colors <- c("#7A8B8B", "#B2DFEE")

ring_width <- as.numeric(0.08)  # To keep the same thickness across all rings
dista <- as.numeric(2)

plot0 <- gheatmap(
  p, data = cc, offset = dista, width = ring_width,
  colnames_angle = 90, font.size = 5
) +
  scale_fill_manual(
    values = cc_colors, name = "Subspecies",
    guide = guide_legend(order = 1)
  )

# Ring 2: Serovar
cc <- data.frame(as.character(meta_filt$serovar), row.names = tree$tip.label)
cc$as.character.meta_filt.serovar. <- gsub("II", "", cc$as.character.meta_filt.serovar.)
colnames(cc) <- "Serovar"
cc_colors <- c("#bbffc5", "#5e8063")

distance <- as.numeric(1)
# Take the latest 'dista' and add 'distance' to push rings away from the center
dista <- distance + dista

plot0 <- plot0 + ggnewscale::new_scale_fill()

plot1 <- gheatmap(
  plot0, data = cc, offset = dista, width = ring_width,
  colnames_angle = 90, font.size = 5
) +
  scale_fill_manual(
    values = cc_colors, name = "Serovar",
    guide = guide_legend(order = 2)
  )

# Ring 3: MLST
cc <- data.frame(as.character(meta_filt$mlst), row.names = tree$tip.label)
colnames(cc) <- "MLST"
cc_colors <- c("#5a5b7d", "#b4b5f9")

distance <- as.numeric(1)
# Take the latest 'dista' and add 'distance' to push rings away from the center
dista <- distance + dista

plot1 <- plot1 + ggnewscale::new_scale_fill()

plot2 <- gheatmap(
  plot1, data = cc, offset = dista, width = ring_width,
  colnames_angle = 90, font.size = 5
) +
  scale_fill_manual(
    values = cc_colors, name = "MLST",
    guide = guide_legend(order = 3)
  )

# Ring 4: Sampling Point
cc <- data.frame(as.character(meta_filt$sample_point), row.names = tree$tip.label)
colnames(cc) <- "Sampling Point"
cc_colors <- c("#00688B", "#FB7185", "#FFDEAD")

distance <- as.numeric(1)
# Take the latest 'dista' and add 'distance' to push rings away from the center
dista <- distance + dista

plot2 <- plot2 + ggnewscale::new_scale_fill()

plot3 <- gheatmap(
  plot2, data = cc, offset = dista, width = ring_width,
  colnames_angle = 90, font.size = 5
) +
  scale_fill_manual(
    values = cc_colors, name = "Sampling Point",
    guide = guide_legend(order = 4)
  ) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13),
    legend.key.size = unit(0.9, "cm"),
    plot.title = element_text(
      hjust = 0.5,  # 0 = left, 0.5 = center, 1 = right
      vjust = 1,    # vertical position (1 is higher)
      size = 14,
      face = "bold"
    )
  )

plot3

# Save fan tree plot
ggsave("./images/cgSNPs_RaxML_kSNP4.pdf", plot = plot3, width = 15, height = 10, dpi = 300)

# High-quality PNG
ggsave(
  "./images/cgSNPs_RaxML_kSNP4.png",
  plot = plot3, width = 15, height = 10, units = "in",
  dpi = 600, device = ragg::agg_png, bg = "white"
)

# TIFF (typical for journals), with LZW compression
ggsave(
  "./images/cgSNPs_RaxML_kSNP4.tiff",
  plot = plot3, width = 15, height = 10, units = "in",
  dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white"
)

# ----------------------------
# hqSNPs trees
# ----------------------------

# hqSNPs
st1815_tree <- read.tree("cfsan/ST1815_ref/RAxML_bestTree.st1815.tree")

st1815 <- ggtree(st1815_tree) +
  geom_tiplab() +
  ggtitle("hqSNPs ST1815") +
  theme_tree2() +
  theme(
    plot.title = element_text(
      hjust = 0.5,  # 0 = left, 0.5 = center, 1 = right
      vjust = 1,    # vertical position (1 is higher)
      size = 14,
      face = "bold"
    )
  )

st1815 <- st1815 +
  xlim(NA, 17.3) +
  geom_tiplab() +
  geom_strip(
    "WR1", "NS4", barsize = 1, color = "#007BFF",
    label = "2-4 hqSNPs", offset.text = .1, offset = 0.85
  ) +
  geom_strip(
    "WR2", "NS3", barsize = 1, color = "#6C7A89",
    label = "0 hqSNPs", offset.text = .1, offset = 0.85
  )

st1815

st3614_tree <- read.tree("cfsan/ST3614_results/RAxML_bestTree.st3614.tree")

st3614 <- ggtree(st3614_tree) +
  geom_tiplab() +
  ggtitle("hqSNPs ST3614") +
  theme_tree2() +
  theme(
    plot.title = element_text(
      hjust = 0.5,  # 0 = left, 0.5 = center, 1 = right
      vjust = 1,    # vertical position (1 is higher)
      size = 14,
      face = "bold"
    )
  )

st3614 <- st3614 +
  xlim(NA, 3.4) +
  geom_tiplab() +
  geom_strip(
    taxa1 = "PBH6", taxa2 = "PBH3",
    barsize = 1, color = "#6C7A89",
    label = "0 hqSNPs",
    offset.text = 0.07, offset = 0.2
  ) +
  geom_strip(
    taxa1 = "PBH5", taxa2 = "PBH1",
    barsize = 1, color = "#007BFF",
    label = "1-2 hqSNPs",
    offset.text = 0.07, offset = 0.2
  )

st3614

# ----------------------------
# Save trees separately
# ----------------------------

# Save tree plots separately in TIFF and PDF
ggsave("cfsan/ST1815_ref/RAxML_bestTree.snp_st1815_tree.pdf", plot = st1815, width = 10, height = 10)
ggsave("cfsan/ST3614_results/RAxML_bestTree.st3614.tree.pdf", plot = st3614, width = 10, height = 10)

# NOTE: The original Rmd saves plot4 here, but plot4 is not defined in the provided code.
# If you intended plot4 to be the core SNP fan tree, use plot3 instead, or define plot4 explicitly.
# Example fix:
# plot4 <- plot3

# If plot4 exists in your full script, keep these lines as-is:
ggsave("kSNP4/coreSNP_kSNP4.tif", plot = plot4, width = 10, height = 10, dpi = 1200)
ggsave("kSNP4/coreSNP_kSNP4.pdf", plot = plot4, width = 10, height = 10)

# ----------------------------
# Combine trees into a single figure
# ----------------------------

# Combine the two hqSNPs trees in a single panel (vertical)
hqSNPs_plots <- ggarrange(
  st1815, st3614,
  ncol = 1, nrow = 2,
  labels = c("B", "C")
)

hqSNPs_plots

# Combine the core tree and the hqSNPs panel
trees_plot <- ggarrange(
  plot4, hqSNPs_plots,
  ncol = 2, nrow = 1,
  labels = c("A", ""),
  widths = c(1.3, 0.7)
)

trees_plot

ggsave("images/all_trees_plot.tif", plot = trees_plot, width = 20, height = 10, dpi = 1200)
ggsave("images/all_trees_plot.pdf", plot = trees_plot, width = 20, height = 10, dpi = 1200)
ggsave("images/all_trees_plot.png", plot = trees_plot, width = 20, height = 10, dpi = 1200)

# ----------------------------
# QUAST table processing
# ----------------------------

# Load QUAST TSV
quast_data <- read.table(
  "assemblies/shovill/assemblies/quast_result/report.tsv",
  header = TRUE, sep = "\t"
)

# Transpose the dataframe
quast_data <- t(quast_data)

# Move the first row to column names
colnames(quast_data) <- quast_data[1, ]

# Remove the first row
quast_data <- quast_data[-1, ]

# Create a column with the Q value
quast_data <- as.data.frame(quast_data)
quast_data$contigs <- as.numeric(quast_data$contigs)
quast_data$`Total length` <- as.numeric(quast_data$`Total length`)

# Q metric: Total length minus (contigs * 1000)
quast_data$Q <- quast_data$`Total length` - (quast_data$contigs * 1000)

# Sort the dataframe by Q (descending)
quast_data <- quast_data[order(quast_data$Q, decreasing = TRUE), ]

# Drop column 22 (as in original script)
quast_data <- quast_data[, -22]
