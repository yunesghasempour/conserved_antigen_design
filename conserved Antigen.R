###############################################
# HPV16 L1 Protein Analysis Pipeline in R
# 
# This script performs the following steps:
# 1. Install and load required packages.
# 2. Query NCBI for HPV16 L1 protein sequences (excluding partials).
# 3. Fetch the sequences in FASTA format and write to file.
# 4. Read sequences and perform multiple sequence alignment (MSA) using ClustalW.
# 5. Compute the consensus sequence from the alignment.
# 6. Calculate a genetic distance matrix and plot it as a heatmap.
# 7. Create segmented MSA plots (in windows of 50 residues) and save them.
# 8. Compute position-wise Shannon entropy from the alignment.
# 9. Plot a phylogenetic tree based on the alignment.
# 10. Generate sequence logo plots for windows (20 amino acids per window) and save them.
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

# Load libraries
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
# 2. Query and Fetch HPV16 L1 Sequences from NCBI
# ------------------------------
# Define query (exclude partial sequences)
query <- "HPV16[Organism] AND L1[Protein Name] NOT partial[Title]"

# Search the 'protein' database for up to 20 results
search_res <- entrez_search(db = "protein", term = query, retmax = 20)

# Fetch the sequences in FASTA format
fasta_data <- entrez_fetch(db = "protein", id = search_res$ids, rettype = "fasta")

# Write the fetched FASTA data to a file
write(fasta_data, file = "HPV16_L1_sequences.fasta")

# ------------------------------
# 3. Read and Align Sequences
# ------------------------------
# Read the FASTA file into an AAStringSet object
seqs <- readAAStringSet("HPV16_L1_sequences.fasta")

# Perform multiple sequence alignment using ClustalW
alignment <- msa(seqs, method = "ClustalW")

# Convert the alignment into a Biostrings AAStringSet
align_biostrings <- as(alignment, "AAStringSet")

# ------------------------------
# 4. Compute Consensus Sequence
# ------------------------------
# Compute the consensus sequence with a threshold of 0.5
consensus_seq <- consensusString(align_biostrings, threshold = 0.5)
cat("\nConsensus Sequence:\n", consensus_seq, "\n")

# Write the consensus sequence to a FASTA file
writeXStringSet(AAStringSet(consensus_seq), "consensus_sequence.fasta")

# ------------------------------
# 5. Genetic Distance Heatmap
# ------------------------------
# Convert the alignment to seqinr alignment object and compute distance matrix
dist_matrix <- dist.alignment(msaConvert(alignment, type = "seqinr::alignment"))

# Convert the distance matrix to a data frame using melt()
df <- melt(as.matrix(dist_matrix))

# Plot the heatmap with ggplot2
ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "thistle") +
  scale_fill_gradient(low = "coral1", high = "pink3") +
  theme_minimal() +
  labs(title = "Genetic Distance Heatmap", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 10, 10, 10))

# ------------------------------
# 6. Create MSA Segment Plots (Window = 50 residues)
# ------------------------------
# Function to create MSA plots for segments of a given window size
create_msa_plots <- function(alignment, window_size = 50) {
  # Get the length of the sequences (assuming all are equal length)
  seq_length <- width(alignment)[1]
  
  # Initialize an empty list to store plots
  plots <- list()
  start_positions <- seq(1, seq_length, by = window_size)
  
  for (i in seq_along(start_positions)) {
    start_pos <- start_positions[i]
    end_pos <- min(start_pos + window_size - 1, seq_length)
    
    # Create the MSA plot for this segment using ggmsa
    p <- ggmsa(alignment, 
               start = start_pos, 
               end = end_pos, 
               color = "Chemistry_AA") +
      ggtitle(sprintf("Positions %d-%d", start_pos, end_pos)) +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    plots[[i]] <- p
  }
  
  return(plots)
}

# Function to save the MSA plots into a PDF and as individual PNG files
save_msa_plots <- function(alignment, filename = "msa_plots.pdf", window_size = 50) {
  plots <- create_msa_plots(alignment, window_size)
  
  # Save all plots in one PDF
  pdf(filename, height = 8, width = 12)
  for (p in plots) {
    print(p)
  }
  dev.off()
  
  # Also save each segment as a separate PNG file in a folder
  dir.create("msa_plots", showWarnings = FALSE)
  for (i in seq_along(plots)) {
    ggsave(
      filename = sprintf("msa_plots/segment_%03d.png", i),
      plot = plots[[i]],
      width = 12,
      height = 6,
      dpi = 300
    )
  }
}

# Save the MSA segment plots (window size = 50 residues)
save_msa_plots(align_biostrings, filename = "HPV16_L1_msa.pdf", window_size = 50)

# Display the MSA plots in R
msa_plots <- create_msa_plots(align_biostrings, window_size = 50)
for (p in msa_plots) {
  print(p)
}

# ------------------------------
# 7. Shannon Entropy Plot for Each Alignment Position
# ------------------------------
# Function to calculate Shannon entropy for a given column of the alignment
calc_entropy <- function(column) {
  # Remove gap characters
  column <- column[column != "-"]
  if (length(column) == 0) return(0)
  freqs <- table(column) / length(column)
  -sum(freqs * log2(freqs + 1e-10))
}

# Convert alignment to a character matrix
align_matrix <- as.matrix(alignment)

# Calculate entropy for each column (position)
entropy_scores <- apply(align_matrix, 2, calc_entropy)

# Plot the entropy scores
plot(entropy_scores, type = "b", pch = 19, col = "blue",
     xlab = "Position in Alignment", ylab = "Shannon Entropy",
     main = "Entropy per Position in MSA")
abline(h = 1, col = "red", lty = 2)  # Threshold line for high conservation

# ------------------------------
# 8. Phylogenetic Tree Plot
# ------------------------------
# Re-calculate distance matrix (if needed) and compute neighbor-joining tree
d <- dist.alignment(msaConvert(alignment, type = "seqinr::alignment"))
tree <- nj(d)
plot(tree, main = "Phylogenetic Tree of HPV16 L1 Sequences")

# ------------------------------
# 9. Create and Save Sequence Logo Plots (Window = 20 residues)
# ------------------------------
# Function to create sequence logo plots for segments of given window size
create_sequence_logos <- function(sequences, window_size = 20, overlap = 0) {
  # If sequences is an AAStringSet, convert to character vector
  if (class(sequences)[1] == "AAStringSet") {
    sequences <- as.character(sequences)
  }
  
  # Get sequence length (assume all sequences have the same length)
  seq_length <- nchar(sequences[1])
  
  # Initialize list to store logo plots
  plots <- list()
  start_positions <- seq(1, seq_length, by = (window_size - overlap))
  
  for (i in seq_along(start_positions)) {
    start_pos <- start_positions[i]
    end_pos <- min(start_pos + window_size - 1, seq_length)
    
    # Extract the subsequences for this window
    subsequences <- sapply(sequences, function(seq) {
      substr(seq, start_pos, end_pos)
    })
    
    # Create sequence logo using ggseqlogo with probability method
    p <- ggseqlogo(subsequences, method = "prob", col_scheme = "chemistry") +
      ggtitle(sprintf("Positions %d-%d", start_pos, end_pos)) +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 11)
      ) +
      scale_y_continuous(limits = c(0, 2))  # Set consistent y-axis limits
    
    plots[[i]] <- p
  }
  
  return(plots)
}

# Function to save the sequence logo plots as a combined PDF and individual PNGs
save_sequence_logos <- function(sequences, filename = "sequence_logos.pdf", 
                                window_size = 20, overlap = 0, ncol = 3) {
  plots <- create_sequence_logos(sequences, window_size, overlap)
  
  # Determine number of rows based on number of plots and desired columns
  n_plots <- length(plots)
  n_rows <- ceiling(n_plots / ncol)
  
  # Save all logos in one PDF
  pdf(filename, height = min(n_rows * 2.5, 20), width = ncol * 6)
  do.call(grid.arrange, c(plots, ncol = ncol))
  dev.off()
  
  # Save each individual logo plot as a PNG file
  dir.create("sequence_logos", showWarnings = FALSE)
  for (i in seq_along(plots)) {
    ggsave(
      filename = sprintf("sequence_logos/window_%03d.png", i),
      plot = plots[[i]],
      width = 6,
      height = 4,
      dpi = 300
    )
  }
}

# Save sequence logo plots for segments of 20 residues (no overlap)
save_sequence_logos(align_biostrings, 
                    filename = "HPV16_L1_logos.pdf",
                    window_size = 20,
                    overlap = 0,
                    ncol = 3)
