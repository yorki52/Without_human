# Without_human
This script automates the analysis of metagenomic data obtained from sequencing samples (cutting human genome info)


The main steps include:

1. **Environment setup**:  
   Before running, it is necessary to download and install the required libraries and tools — Kraken2, KrakenTools, SPAdes, Bracken — as well as download the Kraken database (for example, the standard k2 database, approximately 8 GB). These steps are commented out at the beginning of the notebook for convenience.

2. **Downloading and preparing the Kraken database**:  
   The script automatically downloads and unpacks the Kraken database used for sequence classification.

3. **Processing input data**:  
   In the specified directory, files with paired-end reads (fastq.gz) are located. The script searches for all pairs of R1 and R2 files, checks their existence, and runs Kraken2 on each pair with the database specified. The classification results are saved into individual reports and output files.

4. **Removing human genomic sequences**:  
   During data preparation, human DNA is typically removed to exclude it from analysis since it can occupy a significant portion of sequences and interfere with microbial detection. This is done using specialized tools (e.g., `bbduk`, `bowtie2`, or `kraken2` with a human database). In your case, this step is assumed to be performed beforehand or can be automated.

5. **Combining Kraken reports**:  
   After processing all paired files, the script combines individual Kraken reports into a single consolidated report (`COMBINED.KREPORT`) using the `combine_kreports.py` script. This provides an overview of all sample classifications.

6. **Taxonomic level analysis**:  
   Using the combined report, Bracken estimates the distribution of taxa at various levels — from broad groups (class, order) to more specific ones (family). For each level, statistics are provided: number of taxa with reads above a threshold, below threshold, total reads assigned, and unclassified or unknown reads.

7. **Interpreting results**:  
   The final reports show the microbial composition in the sample as well as the amount of human DNA removed — which is important for assessing sample purity and focusing on microbial components.
