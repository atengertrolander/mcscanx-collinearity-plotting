############################################################
# Convert MCScanX .collinearity output into synteny link file
#
# Input:
#   MCScanX .collinearity file
#
# Output:
#   CSV file with four columns:
#     Chromosome1, Gene1, Chromosome2, Gene2
#
# This output can be used as input for circular synteny plotting
# with circlize.
############################################################


# -----------------------------
# Define input and output files
# -----------------------------
# These are relative paths, so the script assumes you run it from
# the project directory.

input_dir <- "data"
output_dir <- "results"

collinearity_file <- file.path(input_dir, "Cal_Dmel.collinearity.9")
output_link_file <- file.path(output_dir, "link_data.csv")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# -----------------------------
# Read MCScanX collinearity file
# -----------------------------

collinearity_lines <- readLines(collinearity_file)


# -----------------------------
# Initialize an empty list to store gene pairs
# -----------------------------
# Each element of this list will become one row in the final table.

link_list <- list()

# These variables store the chromosome/scaffold names from the current
# alignment block. They are updated each time a new "## Alignment" line
# is encountered.

current_chr1 <- NA
current_chr2 <- NA


# -----------------------------
# Parse the collinearity file
# -----------------------------

for (line in collinearity_lines) {
  
  # Remove leading and trailing whitespace
  line <- trimws(line)
  
  # Skip empty lines
  if (line == "") {
    next
  }
  
  # -----------------------------
  # Identify alignment header lines
  # -----------------------------
  # Example header:
  # ## Alignment 0: score=700.0 e_value=6.6e-35 N=14 NC_004354.4&scaffold_1 plus
  #
  # The chromosome/scaffold pair is the second-to-last field:
  # NC_004354.4&scaffold_1
  
  if (grepl("^## Alignment", line)) {
    
    fields <- strsplit(line, "\\s+")[[1]]
    
    # In standard MCScanX output, the chromosome pair is immediately
    # before the orientation field, for example "plus" or "minus".
    chr_pair <- fields[length(fields) - 1]
    
    # Split chromosome pair into Chromosome1 and Chromosome2
    chr_split <- strsplit(chr_pair, "&")[[1]]
    
    current_chr1 <- chr_split[1]
    current_chr2 <- chr_split[2]
    
    next
  }
  
  
  # -----------------------------
  # Identify gene-pair lines
  # -----------------------------
  # Example gene-pair line:
  # 0-  0:    NP_511100.1    evm.model.scaffold_1.2560    1e-13
  #
  # After splitting on whitespace, the useful fields are:
  #   field 2: Gene1
  #   field 3: Gene2
  #   field 4: e-value, not used here
  
  if (grepl("^[0-9]+-\\s*[0-9]+:", line)) {
    
    fields <- strsplit(line, "\\s+")[[1]]
    
    gene1 <- fields[2]
    gene2 <- fields[3]
    
    link_list[[length(link_list) + 1]] <- data.frame(
      Chromosome1 = current_chr1,
      Gene1 = gene1,
      Chromosome2 = current_chr2,
      Gene2 = gene2,
      stringsAsFactors = FALSE
    )
  }
}


# -----------------------------
# Combine all parsed gene pairs into one data frame
# -----------------------------

link_data <- do.call(rbind, link_list)


# -----------------------------
# Write output link file
# -----------------------------

write.csv(
  link_data,
  file = output_link_file,
  row.names = FALSE,
  quote = FALSE
)

message("Link data written to: ", output_link_file)
message("Number of syntenic gene pairs parsed: ", nrow(link_data))