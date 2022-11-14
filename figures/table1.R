# generating methods table 1
library(ape)
#dat <- read.csv('AM_SAMPLE_AND_DECADE.csv', header=T)
#loc <- read.csv('AMsamples_decimal_coords.csv', header=T)
#dat <- merge(dat, loc, by='Sample.ID')
#write.csv(dat, file='table1.csv', row.names=F)

dat <- read.csv('table1.csv', header=T)

# table S2 for COI dataset
coi <- read.dna('COI653_65pc.fasta', format='fasta', as.character=T)
# 649 samples
#write.csv(rownames(coi),file='tableS2_COI653.csv', row.names=F) modified in excel

CA_coi <- rownames(coi)[grep('CA0',rownames(coi))]


pc5 <- read.dna('5_percent_all_mitogenomes.gatk.nOG.aligned_sname 2.fasta', format='fasta', as.character=T
# old table 1


# new table one from 5pc fasta with CA codes

pc5 <- read.dna('5_percent_all_mitogenomes.gatk.nOG.aligned_F_CODES_forR 2.fasta', format='fasta', as.character=T)
cov <- vector()
for (i in 1:length(pc5[,1])) {
	cov[i] <- 100*(length(which(pc5[i,] %in% c('a', 't', 'g', 'c')))/dim(pc5)[2])
}
tabnew <- cbind(rownames(pc5), cov)
colnames(tabnew) <- c("Sample.ID", 'Percent Coverage')

#write.csv(tabnew, file='table1new.csv', row.names=F)
tabnew <- read.csv(file='table1new.csv', header=T)

tab1_combined <- merge(tabnew, dat, by='Sample.ID', all=T)
write.csv(tab1_combined, file='tab1_combined.csv', row.names=F)

# adding alignment membership
tab1 <- read.csv(file='tab1_combined.csv', header=T)

pc5 <- read.dna('5_percent_all_mitogenomes.gatk.nOG.aligned_F_CODES_forR 2.fasta', format='fasta', as.character=T)
pc65 <- read.dna('65_percent_all_mitogenomes.gatk.nOG.aligned_CODES_sname 2.fasta', format='fasta', as.character=T)
coi <- read.dna('COI653_65pc.fasta', format='fasta', as.character=T)

for (i in 1:length(tab1$header)){
	print(any(grepl(tab1$header[i], rownames(coi))))
	if (any(grepl(tab1$header[i], rownames(coi)))){
		tab1$COI.dataset[i] <- 'Y'
		} else {
		tab1$COI.dataset[i] <- 'N'
		}
}

write.csv(tab1, file='table1_final.csv', row.names=F)