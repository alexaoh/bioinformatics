# Krijnen: Applied Statistics for Bioinformatics using R

library(seqinr)
choosebank()

choosebank("genbank")

# How many ccnd sequences has genbank?
query("ccnd","k=ccnd",virtual=TRUE)$nelem


# How many sequences ccnd sequences has genbank for the species homo sapiens.
query("ccnd3hs","sp=homo sapiens AND k=ccnd3",virtual=TRUE)$nelem


# Let’s download sequences related to the species homo sapiens and a gene name like ”CCND3”.

ccnd3hs <- query("ccnd3hs","sp=homo sapiens AND k=ccnd3@")
ccnd3hs$nelem

sapply(ccnd3hs$req, getName)

sapply(ccnd3hs$req, getLength)

getSequence(ccnd3hs$req[[1]])[1:15]

getAnnot(ccnd3hs$req[[1]])
