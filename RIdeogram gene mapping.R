######################karyotype######################
# Load necessary libraries
library(readxl)
library(tidyverse)
library(RIdeogram)
library(dplyr)
library(readr)
# Read the GFF3 file while skipping metadata lines
gff_data <- read_delim("/Users/macbook/Downloads/Synteny/Hu.gff3", delim = "\t", col_names = FALSE, comment = "#")

# View the structure of the GFF3 data
head(gff_data)

# Split the data into separate columns
colnames(gff_data) <- c("Chr", "Source", "Feature", "Start", "End", "Score", "Strand", "Phase", "Attributes")

# Convert Start and End to numeric
gff_data <- gff_data %>%
  mutate(Start = as.numeric(Start),
         End = as.numeric(End))

# Calculate chromosome lengths based on the max end position of gene features
chrom_lengths <- gff_data %>%
  filter(Feature == "gene") %>%          # Filter for gene entries
  group_by(Chr) %>%                     # Group by chromosome
  summarise(Length = max(End, na.rm = TRUE)) %>% # Calculate length as max end position
  ungroup()

# View the chromosome lengths
print(chrom_lengths)

# Save the results to a CSV or any format
write.csv(chrom_lengths, "karyotype.csv", row.names = FALSE)



######################gene density######################
# Read the gene positions extracted from the GFF3 file
genes <- read.table("/Users/macbook/Downloads/Gene density/genes_positions.txt", header = FALSE, col.names = c("chr", "start", "end"))

# Define the bin size (e.g., 1,000,000 bp = 1 Mb)
bin_size <- 1e6

# Create bins for each chromosome
bins <- genes %>%
  group_by(chr) %>%
  summarize(chr_length = max(end)) %>%
  rowwise() %>%
  do(data.frame(chr = .$chr, 
                start = seq(1, .$chr_length, by = bin_size), 
                end = pmin(seq(bin_size, .$chr_length + bin_size, by = bin_size) - 1, .$chr_length)))

# Count the number of genes in each bin
gene_density <- bins %>%
  rowwise() %>%
  mutate(Value = sum(genes$start >= start & genes$start <= end & genes$chr == chr))

# View the first few rows of the result
head(gene_density)
write.csv(gene_density, "/Users/macbook/Downloads/Gene density/gene_density.csv", row.names = FALSE)

######################gene density and mapping######################
karyotype <- read_excel("/Users/macbook/Downloads/Gene density/karyotype.xlsx")
gene_density <- read_excel("/Users/macbook/Downloads/Gene density/gene_density.xlsx")
wrky_gene <- read_excel("/Users/macbook/Downloads/Gene density/wrky_gene.xlsx")
head(Sme_karyotype)
head(gene_density)
head(wrky_gene)
head(Sme_karyotype)
head(gene_density)
ideogram(karyotype = Sme_karyotype)
convertSVG("chromosome.svg", device = "png")
ideogram(karyotype = Sme_karyotype, overlaid = gene_density)
convertSVG("chromosome.svg", device = "png")
ideogram(karyotype = Sme_karyotype, overlaid = gene_density, label = wrky_gene, label_type = "marker")
convertSVG("chromosome.svg", device = "png")
