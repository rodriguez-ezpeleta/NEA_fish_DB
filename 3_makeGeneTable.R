## REQUIRED LIBRARIES
library(readr)
library(stringr)
library(tidyverse)
library(ape)
library(Orcs)

#Set working directory
setwd("~/OneDrive - AZTI/DB-Cristina/scripts")

data <- read_tsv("output_files/sp_genes.txt", col_names = T)
data <- unique(data)
# Remove UNVERIFIED sequences (https://www.ncbi.nlm.nih.gov/genbank/unverified/)
data <- dplyr::filter(data, !grepl("UNVERIFIED", Sequence_Title))

# Remove erroneous assignations
nrow(data) - length(unique(data$Acc_Number)) #if this value >0 run the step. If =0, move to next step 
n_occur <- data.frame(table(b$Acc_Number))
n_occur<- n_occur[n_occur$Freq > 1, ]
duplicated <- data[data$Acc_Number %in% n_occur$Var1, ]
duplicated$check <- NA
for (i in  c(1:dim(duplicated[1]))) {
  duplicated$check[i]<- ifelse((str_detect(duplicated$Sequence_Title[i],duplicated$Species_Name[i]))== "FALSE", "ERROR", NA)
}
duplicated <- duplicated %>% filter(check =="ERROR")
duplicated$check <-NULL
data <- anti_join(data, duplicated) 

# Assign content of the sequences by descrition (Sequence_Title)
data$regionCOI <- NA
data$regionCOI <- ifelse((str_detect(data$Sequence_Title , "cytocrome c oxidase subunit 1")) == "TRUE"
                          |(str_detect(data$Sequence_Title, "cytochrome oxidase subunit 1")) == "TRUE" 
                          |(str_detect(data$Sequence_Title, "cytochrome c oxidase subunit 1")) == "TRUE" 
                          |(str_detect(data$Sequence_Title,"cytochrome c oxidase 1")) =="TRUE"
                          |(str_detect(data$Sequence_Title,"cytochrome oxidase subunit I gene")) == "TRUE" 
                          |(str_detect(data$Sequence_Title,"cytochrome oxidase subunit I-like gene"))=="TRUE"
                          |(str_detect(data$Sequence_Title, "cytochrome c oxidase subunit I gene")) == "TRUE"
                          |(str_detect(data$Sequence_Title, "cytochrome oxidase subunit I,")) == "TRUE"
                          |(str_detect(data$Sequence_Title, "cytochrome c oxidase subunit I mRNA,")) == "TRUE"
                          |(str_detect(data$Sequence_Title, "cytochrome c oxidase subunits I and")) == "TRUE"
                          |(str_detect(data$Sequence_Title, "cytochrome c oxidase subunit I,")) == "TRUE"
                          |(str_detect(data$Sequence_Title,"cytochrome oxidase I gene")) == "TRUE"
                          |(str_detect(data$Sequence_Title,"COI gene")) =="TRUE"
                          |(str_detect(data$Sequence_Title,"COX1")) == "TRUE"
                          |(str_detect(data$Sequence_Title, "COI (partial)"))== "TRUE"
                          |(str_detect(data$Sequence_Title,"coi gene" ))== "TRUE"
                          |(str_detect(data$Sequence_Title, "cox1 gene"))== "TRUE"
                          |(str_detect(data$Sequence_Title,"(co1)" ))== "TRUE"
                          |(str_detect(data$Sequence_Title, "(COI)"))== "TRUE"
                          |(str_detect(data$Sequence_Title,"(CO1)" ))== "TRUE"
                          |(str_detect(data$Sequence_Title,"(COXI)" ))== "TRUE"
                          |(str_detect(data$Sequence_Title,"(cox1)" ))== "TRUE"
                          |(str_detect(data$Sequence_Title,"(CO-I)" ))== "TRUE"
                          |(str_detect(data$Sequence_Title,"(COX I)" ))== "TRUE"
                          |(str_detect(data$Sequence_Title,"(COX-1)" ))== "TRUE"
                          |(str_detect(data$Sequence_Title,"(coxI)"))== "TRUE"
                          |(str_detect(data$Sequence_Title,"(CoxI)"))== "TRUE"
                          |(str_detect(data$Sequence_Title,"(cox I)" ))== "TRUE",
                          "COI", data$regionCOI)

data$region12S <- NA
data$region12S <- ifelse((str_detect(data$Sequence_Title, "SSU"))== "TRUE" 
                          | (str_detect(data$Sequence_Title, "small subunit ribosomal RNA"))== "TRUE" 
                          | (str_detect(data$Sequence_Title, "12 ribosomal RNA gene"))== "TRUE"
                          | (str_detect(data$Sequence_Title, "12S ribosomal RNA"))== "TRUE"
                          | (str_detect(data$Sequence_Title, "12S rRNA"))== "TRUE"
                          | (str_detect(data$Sequence_Title, "12S"))== "TRUE",
                          "12S", data$region12S)

data$region16S <- NA
data$region16S <- ifelse((str_detect(data$Sequence_Title, "16S ribosomal RNA"))== "TRUE"
                          |(str_detect(data$Sequence_Title, "16S rRNA"))== "TRUE" 
                          |(str_detect(data$Sequence_Title, "16S")) == "TRUE"
                          |(str_detect(data$Sequence_Title, "LSU rRNA"))== "TRUE"
                          |(str_detect(data$Sequence_Title, "large subunit ribosomal RNA"))== "TRUE",
                          "16S", data$region16S) 

data$regioncytb <- NA
data$regioncytb <- ifelse((str_detect(data$Sequence_Title, "cyt b"))== "TRUE" 
                           | (str_detect(data$Sequence_Title, "cytb"))== "TRUE" 
                           | (str_detect(data$Sequence_Title,"CYTB"))== "TRUE"
                           | (str_detect(data$Sequence_Title, "cytochrome B"))== "TRUE"
                           | (str_detect(data$Sequence_Title, "cytochrome b"))== "TRUE"
                           | (str_detect(data$Sequence_Title, "cytochrome oxidase B")) == "TRUE"
                           | (str_detect(data$Sequence_Title, "cytochrome-b"))== "TRUE",
                           "cytb", data$regioncytb)

data$wholemtDNA <- NA
data$wholemtDNA <- ifelse((str_detect(data$Sequence_Title, "complete genome"))== "TRUE" 
                          | (str_detect(data$Sequence_Title, "complete mitochondiral genome"))== "TRUE" 
                          | (str_detect(data$Sequence_Title,"complete mitochondrial DNA sequence"))== "TRUE"
                          | (str_detect(data$Sequence_Title, "complete mitochondrial DNA genome"))== "TRUE"
                          | (str_detect(data$Sequence_Title, "genome assembly, organelle: mitochondrion"))== "TRUE"
                          | (str_detect(data$Sequence_Title, "mitochondrial DNA, almost complete sequence")) == "TRUE"
                          | (str_detect(data$Sequence_Title, "mitochondrial DNA, compete genome"))== "TRUE"
                          | (str_detect(data$Sequence_Title, "mitochondrial DNA, complete and partial sequence"))== "TRUE"
                          | (str_detect(data$Sequence_Title, "mitochondrial DNA, complete sequence"))== "TRUE"
                          | (str_detect(data$Sequence_Title, "mitochondrial genome"))== "TRUE"
                          | (str_detect(data$Sequence_Title, "mitochondrion, complete sequence"))== "TRUE"
                          | (str_detect(data$Sequence_Title, "mitochondrionial geneome"))== "TRUE"
                          | (str_detect(data$Sequence_Title, "whole genome shotgun"))== "TRUE",
                          "wholemtDNA", data$wholemtDNA)

for (i in c(1:dim(data)[1])) {
  if (!is.na(data$wholemtDNA[i])) {
    data$regionCOI[i] <- "COI"
    data$region12S[i] <- "12S"
    data$region16S[i] <- "16S"
    data$regioncytb[i] <- "cytb"
  }
}
data$wholemtDNA <- NULL

data$partialmtDNA <- NA
data$partialmtDNA <- ifelse((str_detect(data$Sequence_Title, "mitochondrion, partiale genome"))== "TRUE" 
                            | (str_detect(data$Sequence_Title, "mitochondrial DNA, partial sequence"))== "TRUE"
                            | (str_detect(data$Sequence_Title, "mitochondrial DNA, partial genome"))== "TRUE"
                            | (str_detect(data$Sequence_Title, "mitochondrion, partial genome"))== "TRUE",
                            "partialmtDNA", data$partialmtDNA)

data$othermtDNA <- NA
data$othermtDNA <- ifelse(is.na(data$regionCOI) & is.na(data$region12S) & is.na(data$region16S) & is.na(data$regioncytb) & is.na(data$partialmtDNA),
                          "other mtDNA", data$othermtDNA)

# In order to identify partialmtDNA, we will download the sequences to perform an alignment 
# Create a vector with GenBank accession numbers
partialmtDNA <- data %>% filter(partialmtDNA == "partialmtDNA" )
partialmtDNA_accession_numbers <- unlist(partialmtDNA[,4])
# Download FASTA file
partialmtDNA_sequences <- read.GenBank(partialmtDNA_accession_numbers, quiet = F) 
write.dna(partialmtDNA_sequences, file ="PartialmtDNA.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10) #BLAST input file

#Among the already identified sequences, choose complete sequences corresponding to target regions to be used as references for the alignment 
data$Sequence_Length <- as.numeric(data$Sequence_Length)
gene_COI <- data %>% filter(regionCOI=="COI" & is.na(region12S) & is.na(region16S) & is.na(regioncytb))
gene_COI <- filter(gene_COI, Sequence_Length > 1200 & Sequence_Length < 1800)
gene_COI <- dplyr::filter(gene_COI, !grepl("partial", Sequence_Title))
gene_COI <- gene_COI$Acc_Number
ref_seq_COI <- read.GenBank(gene_COI,quiet = F)
write.dna(ref_seq_COI, file ="completeCOI_ref_seq.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10) #BLAST input file (COI alignment)

gene_12S <- data %>% filter(is.na(regionCOI) & region12S=="12S" & is.na(region16S) & is.na(regioncytb))
gene_12S <- filter(gene_12S, Sequence_Length > 700 & Sequence_Length < 1300) #if & does not work, try |
gene_12S <- dplyr::filter(gene_12S, !grepl("partial", Sequence_Title))
gene_12S <- gene_12S$Acc_Number
ref_seq_12S <- read.GenBank(gene_12S,quiet = F)
write.dna(ref_seq_12S, file ="complete12S_ref_seq.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10) #BLAST input file (12S alignment))

gene_16S <- data %>% filter(is.na(regionCOI) & is.na(region12S) & region16S=="16S" &  is.na(regioncytb))
gene_16S <- filter(gene_16S, Sequence_Length > 1200 & Sequence_Length < 1800) 
gene_16S <- dplyr::filter(gene_16S, !grepl("partial", Sequence_Title))
gene_16S <- gene_16S$Acc_Number
ref_seq_16S <- read.GenBank(gene_16S,quiet = F)
write.dna(ref_seq_16S, file ="complete16S_ref_seq.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10) #BLAST input file (16S alignment)

gene_cytb <- data %>% filter(is.na(regionCOI) & is.na(region12S) & is.na(region16S) & regioncytb=="cytb")
gene_cytb <- dplyr::filter(gene_cytb, Sequence_Length > 800 & Sequence_Length < 1400)
gene_cytb <- dplyr::filter(gene_cytb, !grepl("partial", Sequence_Title))
gene_cytb <- gene_cytb[sample(nrow(gene_cytb), 100), ]#random sampling of 100 sequences
gene_cytb <- gene_cytb$Acc_Number
ref_seq_cytb <- read.GenBank(gene_cytb,quiet = F)
write.dna(ref_seq_cytb, file ="completecytb_ref_seq.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10) #BLAST input file (cytb alignment)

#Make a database with PartialmtDNA fasta sequences
cmd: makeblastdb -in PartialmtDNA.fasta -dbtype nucl -out PartialmtDNA_Database

#For each gene, perform an aligment by Basic Local Aligment Seach Tool with the created FASTA file
cmd: blastn -db PartialmtDNA_Database -query completeCOI_ref_seq.fasta -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue stitle' -perc_identity 60 -evalue 0.0001 -out Output_blast_COI_PartialmtDNA.txt
cmd: blastn -db PartialmtDNA_Database -query complete12S_ref_seq.fasta -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue stitle' -perc_identity 60 -evalue 0.0001 -out Output_blast_12S_PartialmtDNA.txt
cmd:blastn -db PartialmtDNA_Database -query complete16S_ref_seq.fasta -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue stitle' -perc_identity 60 -evalue 0.0001 -out Output_blast_16S_PartialmtDNA.txt
cmd:blastn -db PartialmtDNA_Database -query completecytb_ref_seq.fasta -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue stitle' -perc_identity 60 -evalue 0.0001 -out Output_blast_cytb_PartialmtDNA.txt

blast_COI <- read.delim("Output_blast_COI_PartialmtDNA.txt", header=FALSE)
blast_COI$regionCOI <- c(rep("COI", nrow(blast_COI)))
blast_COI <- dplyr::select(blast_COI, "V2", "regionCOI")
blast_COI <- unique(blast_COI)

blast_12S <- read.delim("Output_blast_12S_PartialmtDNA.txt", header=FALSE)
blast_12S$region12S <- c(rep("12S", nrow(blast_12S)))
blast_12S <- dplyr::select(blast_12S, "V2", "region12S")
blast_12S <- unique(blast_12S)

blast_16S <- read.delim("Output_blast_16S_PartialmtDNA.txt", header=FALSE)
blast_16S$region16S <- c(rep("16S", nrow(blast_16S)))
blast_16S <- dplyr::select(blast_16S, "V2", "region16S")
blast_16S <- unique(blast_16S)

blast_cytb <- read.delim("Output_blast_cytb_PartialmtDNA.txt", header=FALSE)
blast_cytb$regioncytb <- c(rep("cytb", nrow(blast_cytb)))
blast_cytb <- dplyr::select(blast_cytb, "V2", "regioncytb")
blast_cytb <- unique(blast_cytb)

blast <- Orcs::merge(list(blast_COI, blast_12S, blast_16S, blast_cytb), by = "V2", all=T)
names(blast)[1] <- "Acc_Number"

# Add new identified sequences to the main dataset
data <- left_join(data, blast, by = c("Acc_Number")) %>% 
  mutate(regionCOI = ifelse(is.na(regionCOI.x), regionCOI.y, regionCOI.x)) %>% 
  mutate(region12S = ifelse(is.na(region12S.x), region12S.y, region12S.x)) %>% 
  mutate(region16S = ifelse(is.na(region16S.x), region16S.y, region16S.x)) %>% 
  mutate(regioncytb = ifelse(is.na(regioncytb.x), regioncytb.y, regioncytb.x)) %>% 
  dplyr::select(Species_Name,Taxa_Id, Gi, Acc_Number, Sequence_Title, Sequence_Length, regionCOI, region12S, region16S, regioncytb, partialmtDNA, othermtDNA)
data$partialmtDNA <- NULL
# Remove other mtDNA sequences
data <- data %>% filter(is.na(othermtDNA))
data$othermtDNA<-NULL

# Save dataset
write.table(data, "Genbank_available_sequences.txt", row.names = FALSE, col.names = TRUE)

# Download COI and 12S gene sequences (will be used later)
Sequences_COI <- data %>% filter(data$regionCOI=="COI")
Sequences_COI <- unique(Sequences_COI$Acc_Number)
Sequences_COI <- read.GenBank(Sequences_COI, chunk.size = 300, quiet=F) 
write.dna(Sequences_COI, file ="Sequences_COI.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10) 

Sequences_12S <- data %>% filter(data$region12S=="12S")
Sequences_12S <- unique(Sequences_12S$Acc_Number)
Sequences_12S <- read.GenBank(Sequences_12S, chunk.size = 300, quiet=F) 
write.dna(Sequences_12S, file ="Sequences_12S.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10) 

# Make summary file
# Prepare data
data_elongated <- data %>% gather(colname, Gene, -c(Species_Name, Taxa_Id, Gi, Acc_Number, Sequence_Title, Sequence_Length))
data_elongated <- data_elongated %>% filter(!is.na(Gene))
data_elongated <- data_elongated %>% dplyr::select(Species_Name, Gene)

species <- unique(data_elongated$Species_Name)
genes <- unique(data_elongated$Gene)
vector<-c()
for (i in species) { 
  for (ii in genes) {
    vector <-c(vector,length(data_elongated[(data_elongated$Species_Name==i) & (data_elongated$Gene == ii), 2]))
  }
}
summary <- data.frame(matrix(data = vector, byrow=TRUE, nrow=length(species), ncol=length(genes)))
colnames(summary) <- genes
rownames(summary) <- species
summary$Species_Name <- rownames(summary)

# Import taxonomy file created in previous step
Species_list <- read.csv("output_files/Species_list.txt", sep = "")

# Merge taxonomy with gene information
data_with_taxa <- merge(Species_list, summary, by="Species_Name", all = T) 
col_order <- c("Kingdom", "Phylum", "Class_Gigaclass", "Order", "Family", "Genus", "Species_Name", "COI", "12S", "16S", "cytb")
data_with_taxa <- data_with_taxa[, col_order]
data_with_taxa[is.na(data_with_taxa)]<- 0

# Export data
write.table(data_with_taxa, "Available_sequences_marine_fish.txt", sep = " ", dec = ".",
            row.names = FALSE, col.names = TRUE)
