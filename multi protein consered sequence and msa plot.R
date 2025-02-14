###############################################
# Multi-Protein Consensus & MSA Chart Pipeline
#
# For a given strain and protein type, this script:
#   1. Queries NCBI for protein sequences (excluding partials)
#   2. Saves fetched sequences (FASTA) in a folder
#   3. Performs multiple sequence alignment (MSA) using ClustalW
#   4. Computes the consensus sequence and replaces '?' with 'X'
#   5. Saves the consensus sequence in a FASTA file
#   6. Creates segmented MSA plots (50 residues per segment)
#   7. Saves all segmented MSA plots as a PDF and as individual PNGs
#
# All outputs are saved in a folder named "strain_proteinType"
###############################################

# ------------------------------
# 1. Install and Load Required Packages
# ------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("rentrez", quietly = TRUE))
  BiocManager::install("rentrez")
if (!requireNamespace("msa", quietly = TRUE))
  BiocManager::install("msa")
if (!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")
if (!requireNamespace("ggmsa", quietly = TRUE))
  install.packages("ggmsa")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("ape", quietly = TRUE))
  install.packages("ape")
if (!requireNamespace("reshape2", quietly = TRUE))
  install.packages("reshape2")
if (!requireNamespace("gridExtra", quietly = TRUE))
  install.packages("gridExtra")
if (!requireNamespace("ggseqlogo", quietly = TRUE))
  install.packages("ggseqlogo")

library(rentrez)
library(msa)
library(Biostrings)
library(ggmsa)
library(ggplot2)
library(ape)
library(reshape2)
library(gridExtra)
library(ggseqlogo)

# ------------------------------
# 2. Define the Pipeline Function
# ------------------------------
run_protein_pipeline <- function(strain, proteinType, retmax = 20, window_size_msa = 50) {
  # Create output folder based on strain and proteinType (e.g., "HPV16_L1")
  folderName <- sprintf("%s_%s", strain, proteinType)
  if (!dir.exists(folderName)) dir.create(folderName)
  
  # ------------------------------
  # 3. Query NCBI and Fetch Sequences
  # ------------------------------
  query <- sprintf('%s[Organism] AND %s[Protein Name] NOT partial[Title]', strain, proteinType)
  search_res <- entrez_search(db = "protein", term = query, retmax = retmax)
  fasta_data <- entrez_fetch(db = "protein", id = search_res$ids, rettype = "fasta")
  
  # Save FASTA data
  fasta_file <- file.path(folderName, sprintf("%s_%s_sequences.fasta", strain, proteinType))
  write(fasta_data, file = fasta_file)
  
  # ------------------------------
  # 4. Read Sequences and Perform MSA
  # ------------------------------
  seqs <- readAAStringSet(fasta_file)
  alignment <- msa(seqs, method = "ClustalW")
  align_biostrings <- as(alignment, "AAStringSet")
  
  # ------------------------------
  # 5. Compute Consensus Sequence and Fix Invalid Characters
  # ------------------------------
  consensus_seq <- consensusString(align_biostrings, threshold = 0.5)
  # Replace any '?' characters with 'X'
  consensus_seq_fixed <- gsub("\\?", "X", consensus_seq)
  cat(sprintf("\nConsensus for %s %s:\n%s\n", strain, proteinType, consensus_seq_fixed))
  
  consensus_file <- file.path(folderName, sprintf("consensus_%s_%s.fasta", strain, proteinType))
  writeXStringSet(AAStringSet(consensus_seq_fixed), consensus_file)
  
  # ------------------------------
  # 6. Create Segmented MSA Plots (Window = 50 residues)
  # ------------------------------
  create_msa_plots <- function(alignment, window_size) {
    seq_length <- width(alignment)[1]  # assume equal length for all sequences
    plots <- list()
    start_positions <- seq(1, seq_length, by = window_size)
    
    for (i in seq_along(start_positions)) {
      start_pos <- start_positions[i]
      end_pos <- min(start_pos + window_size - 1, seq_length)
      p <- ggmsa(alignment, start = start_pos, end = end_pos, color = "Chemistry_AA") +
        ggtitle(sprintf("Positions %d-%d", start_pos, end_pos)) +
        theme(plot.title = element_text(size = 12, face = "bold"))
      plots[[i]] <- p
    }
    return(plots)
  }
  
  # Save segmented MSA plots as PDF and PNGs
  save_msa_plots <- function(alignment, folder, window_size) {
    plots <- create_msa_plots(alignment, window_size)
    
    # Save all plots in one PDF
    pdf_file <- file.path(folder, sprintf("%s_%s_MSA_segments.pdf", strain, proteinType))
    pdf(pdf_file, height = 8, width = 12)
    for (p in plots) print(p)
    dev.off()
    
    # Save individual segments as PNG
    msa_dir <- file.path(folder, "msa_plots")
    if (!dir.exists(msa_dir)) dir.create(msa_dir)
    for (i in seq_along(plots)) {
      ggsave(filename = file.path(msa_dir, sprintf("segment_%03d.png", i)),
             plot = plots[[i]], width = 12, height = 6, dpi = 300)
    }
  }
  
  save_msa_plots(align_biostrings, folderName, window_size_msa)
  
  cat(sprintf("Pipeline completed for %s %s. All outputs are saved in folder: %s\n", strain, proteinType, folderName))
}


# ------------------------------
# 3. Run the Pipeline for Multiple Proteins
# ------------------------------
# Define a list of proteins with strain and protein type
proteinList <- list(
  list(strain = "HPV16", proteinType = "L2"),
  list(strain = "HPV16", proteinType = "L1"),
  list(strain = "HPV16", proteinType = "E6"),
  list(strain = "HPV16", proteinType = "E7")
)

for (protein in proteinList) {
  run_protein_pipeline(strain = protein$strain, proteinType = protein$proteinType)
}
