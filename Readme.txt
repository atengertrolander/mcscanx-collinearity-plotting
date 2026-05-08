# Synteny analysis and circular plotting scripts

This repository contains input files and custom R scripts used to parse MCScanX collinearity output and generate circular synteny plots with the R package `circlize`.

The synteny analysis itself is described in the associated manuscript Methods. The scripts provided here are intended to document the downstream processing and visualization steps used to convert MCScanX collinearity output into a simplified link table and plot syntenic relationships between chromosomes/scaffolds.

## Repository structure

mcscanx-collinearity-plotting/
├── README.md
├── 01_parsing_mcscanx_collinearity_data.R
├── 02_plot_synteny.R
├── data/
│   ├── Cal_Dmel.collinearity.9
│   ├── Cal_Dmel.gff
│   ├── Scaff_lengths.csv
└── results/
    ├── link_data.csv
    └── synteny_circos_plot.pdf

R scripts:
1. 01_parsing_mcscanx_collinearity_data.R converts the MCScanX .collinearity output into a simplified link table by extracting each alignment block’s chromosome/scaffold pair and the syntenic gene pairs listed within that block.
2. 02_plot_synteny.R uses the parsed synteny link table, gene coordinate file, scaffold length file, and optional highlight regions to generate a circular synteny plot with the R package circlize

Data files:
1. Cal_Dmel.collinearity.9 is the alignment file that McScanX outputs
2. Cal_Dmel.gff is a simplified gene coordinate file used to assign genomic positions to the syntenic gene pairs. Although this file has a `.gff` extension, it is **not a standard nine-column GFF3 file**. It is a simplified four-column gene position file in the format expected by the plotting script and similar to the BED/GFF-like coordinate file used by MCScanX.
3. Scaff_lengths.csv contains the chromosome/scaffold IDs and coordinate ranges used to initialize the circular plot layout, with one row per chromosome or scaffold and columns specifying the sequence name, start coordinate, and end coordinate.

