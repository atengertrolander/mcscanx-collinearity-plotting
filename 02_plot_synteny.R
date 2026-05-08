############################################################
# Plot synteny relationships as a circular plot
#
# This script takes:
#   1. scaffold/chromosome length information
#   2. MCScanX-derived collinearity link data
#   3. a BED/GFF-like file with gene coordinates
#   4. optional genomic regions to highlight
#
# It produces a circular synteny plot using the R package circlize.
############################################################

# -----------------------------
# Load required libraries
# -----------------------------

library(circlize)


# -----------------------------
# Define input and output paths
# -----------------------------
# All paths are relative to the folder from which this script is run.
# To run the script, place this file in the project directory or update
# the file names below to match your directory structure.

input_dir <- "data"
output_dir <- "results"

chr_lengths_file <- file.path(input_dir, "Scaff_lengths.csv")
gene_coordinate_file <- file.path(input_dir, "Cal_Dmel.gff")
link_data_file <- file.path(input_dir, "link_data.csv")
highlight_regions_file <- file.path(input_dir, "Highlight_regions.csv")

output_plot_file <- file.path(output_dir, "synteny_circos_plot.pdf")

# Create the output directory if it does not already exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# -----------------------------
# Load chromosome/scaffold length data
# -----------------------------
# This file should contain one row per chromosome or scaffold.
# Expected columns include:
#   Chromosome: chromosome/scaffold ID
#   Start: start coordinate, usually 0 or 1
#   End: chromosome/scaffold length

chr_data <- read.csv(chr_lengths_file, header = TRUE, stringsAsFactors = FALSE)


# -----------------------------
# Load gene coordinate data
# -----------------------------
# This file should contain gene positions used to place synteny links.
# Expected columns, in order:
#   Chromosome
#   GeneID
#   Start
#   End
#
# Although the file name here ends in .gff, this script expects a
# simplified four-column table rather than a standard nine-column GFF3.

bed_file <- read.table(
  gene_coordinate_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(bed_file) <- c("Chromosome", "GeneID", "Start", "End")


# -----------------------------
# Load synteny link data
# -----------------------------
# This file should contain pairs of syntenic genes extracted from the
# MCScanX collinearity output.
#
# Expected columns:
#   Chromosome1
#   Gene1
#   Chromosome2
#   Gene2
#
# If the input CSV has different column names, they are standardized below.

link_data <- read.csv(link_data_file, header = TRUE, stringsAsFactors = FALSE)

colnames(link_data) <- c("Chromosome1", "Gene1", "Chromosome2", "Gene2")


# -----------------------------
# Add genomic coordinates to each syntenic gene pair
# -----------------------------
# The synteny link file contains gene IDs, but circos.link() needs
# genomic coordinates. These merges add the start and end positions for
# Gene1 and Gene2 from the gene coordinate table.

link_data <- merge(
  link_data,
  bed_file,
  by.x = "Gene1",
  by.y = "GeneID",
  all.x = TRUE
)

link_data <- merge(
  link_data,
  bed_file,
  by.x = "Gene2",
  by.y = "GeneID",
  suffixes = c(".1", ".2")
)

# Remove duplicate chromosome columns introduced by the merge.
# Chromosome1 and Chromosome2 from the original link file are retained.
link_data <- link_data[, !(names(link_data) %in% c("Chromosome.1", "Chromosome.2"))]


# -----------------------------
# Optional quality check
# -----------------------------
# Check whether any syntenic genes failed to match the coordinate file.
# Missing coordinates will prevent those links from being plotted correctly.

missing_coordinates <- link_data[
  is.na(link_data$Start.1) |
    is.na(link_data$End.1) |
    is.na(link_data$Start.2) |
    is.na(link_data$End.2),
]

if (nrow(missing_coordinates) > 0) {
  warning(
    nrow(missing_coordinates),
    " synteny links have missing coordinates and may not plot correctly."
  )
}


# -----------------------------
# Define colors for synteny links
# -----------------------------
# Links are colored according to Chromosome1.
# Any chromosome not listed here will be plotted in black.

chromosome_colors <- c(
  "NC_004354.4" = "#D36027",
  "NT_033779.5" = "#F0E442",
  "NT_033778.4" = "#009E73",
  "NT_037436.4" = "#0073B2",
  "NT_033777.3" = "#5AB4E5",
  "NC_004353.4" = "#C1C1C1"
)


# -----------------------------
# Open output plotting device
# -----------------------------
# The plot will be written directly to a PDF file.

pdf(output_plot_file, width = 10, height = 10)


# -----------------------------
# Initialize circular plot layout
# -----------------------------
# circos.initialize() creates one sector for each chromosome/scaffold.
# The x-axis limits for each sector are defined by the start and end
# coordinates in chr_data.

circos.par(cell.padding = c(0.02, 0, 0.02, 0))

circos.initialize(
  factors = chr_data$Chromosome,
  xlim = chr_data[, c("Start", "End")]
)


# -----------------------------
# Add chromosome/scaffold labels
# -----------------------------
# This track labels each sector using the chromosome/scaffold ID.

circos.trackPlotRegion(
  factors = chr_data$Chromosome,
  ylim = c(0, 1),
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      0.5,
      CELL_META$sector.index,
      facing = "inside",
      niceFacing = TRUE
    )
  }
)


# -----------------------------
# Highlight selected genomic regions
# -----------------------------
# This optional input file should contain regions to highlight.
# Expected columns:
#   SeqID: chromosome/scaffold ID
#   Start: start coordinate
#   End: end coordinate
#   Color: optional color for the highlighted region
#
# If no Color column is present, highlighted regions are shown in red.

if (file.exists(highlight_regions_file)) {
  highlight_regions <- read.csv(
    highlight_regions_file,
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(nrow(highlight_regions))) {
    highlight_color <- if ("Color" %in% colnames(highlight_regions)) {
      highlight_regions$Color[i]
    } else {
      "red"
    }
    
    circos.rect(
      xleft = highlight_regions$Start[i],
      ybottom = 0,
      xright = highlight_regions$End[i],
      ytop = 1,
      sector.index = highlight_regions$SeqID[i],
      col = highlight_color,
      border = "black"
    )
  }
}


# -----------------------------
# Plot synteny links
# -----------------------------
# Each row of link_data represents one syntenic gene pair.
# circos.link() draws a curved link between the genomic coordinates
# of Gene1 and Gene2.
#
# Links are colored by Chromosome1 using the color map defined above.

for (i in seq_len(nrow(link_data))) {
  
  link_color <- chromosome_colors[link_data$Chromosome1[i]]
  
  if (is.na(link_color)) {
    link_color <- "#000000"
  }
  
  circos.link(
    sector.index1 = link_data$Chromosome1[i],
    point1 = c(link_data$Start.1[i], link_data$End.1[i]),
    sector.index2 = link_data$Chromosome2[i],
    point2 = c(link_data$Start.2[i], link_data$End.2[i]),
    col = link_color
  )
}


# -----------------------------
# Finish plot
# -----------------------------

circos.clear()
dev.off()

message("Synteny circos plot written to: ", output_plot_file)