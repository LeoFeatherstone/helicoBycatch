# script converts tips labs in 65% fasta to match the dates.txt file

library(ape)

aln <- read.dna("65percent_mitogenomes.fasta", format = "fasta")

rownames(aln)
rownames(aln) <- gsub(rownames(aln), pattern = "AM", replacement = "AUS_AM")
rownames(aln)
rownames(aln) <- gsub(rownames(aln), pattern = "Au", replacement = "AUS")
rownames(aln) <- gsub(rownames(aln), pattern = "Br", replacement = "BRA")
rownames(aln) <- gsub(rownames(aln), pattern = "Ch", replacement = "CHN")
rownames(aln) <- gsub(rownames(aln), pattern = "Ug", replacement = "UGA")
rownames(aln) <- gsub(rownames(aln), pattern = "In", replacement = "IND")
rownames(aln) <- gsub(rownames(aln), pattern = "Se", replacement = "SEN")
rownames(aln) <- gsub(rownames(aln), pattern = "Nz", replacement = "NZL")
rownames(aln) <- gsub(rownames(aln), pattern = "Fr", replacement = "FRA")
rownames(aln) <- gsub(rownames(aln), pattern = "Ma", replacement = "MDG")
rownames(aln) <- gsub(rownames(aln), pattern = "_..s", replacement = "")

# check
dates <- read.delim("65%_ages.txt", head = FALSE)
all(rownames(aln) %in% dates[, 1]) & all(dates[, 1] %in% rownames(aln))

write.dna(aln, file = "65percent_mitogenomes.fasta", format = "fasta")
