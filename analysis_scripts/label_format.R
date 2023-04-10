# script to append dates to fasta headers
library(ape)
aln_files <- dir(pattern = ".+.fasta")
aln <- lapply(aln_files, function(x) read.dna(x, format = "fasta"))
aln <- lapply(aln_files, function(x) read.FASTA(x, type = "DNA"))
names(aln) <- aln_files
date_files <- dir(pattern = ".+.txt")
dates <- lapply(date_files, function(x) read.delim(x, header = FALSE))
names(dates) <- aln_files

# reformat .fasta headers
for (i in seq_along(aln)) {
    for (j in seq_along(rownames(aln[[i]]))) {
        which <- grep(dates[[i]][, 1], pattern = rownames(aln[[i]])[j])
        tmp <- dates[[i]][which, 2]
        rownames(aln[[i]])[j] <- paste0(rownames(aln[[i]])[j], "|", tmp)
    }
}
