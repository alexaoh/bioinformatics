# Biostrings functions

# Example:
# Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
# NC_045512 29903 bp ss-RNA     linear 

setwd("/home/ajo/gitRepos/bioinformatics/Task3")
library(Biostrings)

file1 <- "NC_045512.2.fasta"

seqDNA <- readDNAStringSet(file1)  # DNA

# seqDNA features

seqDNA

names(seqDNA)
width(seqDNA)
length(seqDNA) # Numero de secuencias dentro de seqDNA
toString(seqDNA)


# extract several segments of DNA
at1 <- IRanges(start = c(266,13468), end = c(13468,21555))
gene_wuhan <- extractAt(seqDNA,at1)


# DNAStringSet
unlist(unlist(extractAt(seqDNA,at1)))


# Coding Sequence (CDS): Protein

CDS_wuhan <- AAStringSet(translate(unlist(unlist(extractAt(seqDNA,at1)))))


# eliminate stop codon

CDS_wuhan <- subseq(CDS_wuhan, start = 1, end= width(CDS_wuhan)-1)



# orf1a polyprotein [Wuhan seafood market pneumonia virus]
# NCBI Reference Sequence: YP_009724389.1

file2 = "YP_009724389.1.fasta"

seqAA <- readAAStringSet(file2)  # AA

# seqDNA features

seqAA

# compare CDS that it was obtained in the program between CDS in the database 

toString(seqAA) == toString(CDS_wuhan)

