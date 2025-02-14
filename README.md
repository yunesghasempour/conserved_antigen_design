HPV16 L1 Protein Analysis Pipeline
This repository contains an R script that performs a comprehensive analysis of HPV16 L1 protein sequences retrieved from NCBI. The script encompasses the following key steps:

Data Retrieval: Query and fetch HPV16 L1 protein sequences (excluding partial sequences) from NCBI.
Multiple Sequence Alignment (MSA): Align the retrieved sequences using ClustalW.
Consensus Sequence: Compute the consensus sequence from the alignment.
Genetic Distance Heatmap: Calculate a genetic distance matrix and visualize it as a heatmap.
MSA Segmentation: Create and save segmented MSA plots (in windows of 50 residues) to better inspect long alignments.
Entropy Analysis: Calculate and plot the Shannon entropy for each alignment position.
Phylogenetic Analysis: Construct and plot a phylogenetic tree using the neighbor-joining method.
Sequence Logos: Generate sequence logo plots (for windows of 20 residues) using ggseqlogo, and save them both as a combined PDF and as individual PNG files.
Requirements
The script requires the following R packages:

BiocManager
rentrez
msa
Biostrings
ggmsa
ggplot2
ape
reshape2
gridExtra
ggseqlogo
The script includes code to automatically install any missing packages.

Installation
Clone or download the repository, then open the HPV16_L1_Analysis.R script in RStudio or run it from the R console:

bash
Copy
Edit
git clone https://github.com/yourusername/HPV16_L1_Analysis.git
Usage
Simply run the script in R. The script will perform the following actions:

Fetch Data: Retrieve HPV16 L1 protein sequences from NCBI (excluding partials) and save them as HPV16_L1_sequences.fasta.
Alignment: Read the sequences, perform MSA using ClustalW, and convert the result to a Biostrings object.
Consensus: Compute and save the consensus sequence to consensus_sequence.fasta.
Heatmap: Calculate the genetic distance matrix, convert it into a data frame, and plot a heatmap with ggplot2.
MSA Segmentation: Create MSA plots for segments of 50 residues, save all segments in a single PDF (HPV16_L1_msa.pdf), and also export individual PNG files in the msa_plots/ folder.
Entropy Plot: Compute the Shannon entropy for each alignment position and display a plot with a conservation threshold.
Phylogenetic Tree: Generate a neighbor-joining tree from the alignment and plot it.
Sequence Logos: Divide the aligned sequences into windows of 20 residues (with no overlap by default) and generate sequence logos using ggseqlogo. These logos are saved in a PDF (HPV16_L1_logos.pdf) and as individual PNG files in the sequence_logos/ directory.
Output Files
After running the script, you will find the following files in your working directory:

HPV16_L1_sequences.fasta: Fetched sequences in FASTA format.
consensus_sequence.fasta: Consensus sequence derived from the alignment.
HPV16_L1_msa.pdf: PDF containing segmented MSA plots.
msa_plots/: Folder with individual segment PNG plots.
Genetic Distance Heatmap: Displayed in R.
Entropy Plot: Displayed in R.
Phylogenetic Tree: Displayed in R.
HPV16_L1_logos.pdf: PDF with sequence logo plots.
sequence_logos/: Folder with individual sequence logo PNG files.
Customization
Parameters such as window sizes, overlap, color schemes, and output dimensions can be easily adjusted by modifying the corresponding parameters in the script functions:

MSA Segmentation: Change the window_size parameter in the create_msa_plots() function.
Sequence Logos: Modify window_size, overlap, and ncol parameters in the save_sequence_logos() function.
Plot Aesthetics: Adjust ggplot2 themes and scales within the plotting sections to suit your preferences.
License
This project is licensed under the MIT License. See the LICENSE file for details.

Contact
For questions or suggestions, please open an issue in this repository or contact [Your Name] at [your.email@example.com].

You can now commit this README.md along with your R script to your GitHub repository. Enjoy sharing your work!
