library(ggplot2)
library(reshape2)
library(cowplot)
library(gridGraphics)
library(plot.matrix)
library(ape)

setwd(paste0(getwd(), "/figures"))

## data for decay plot
dat <- read.csv("table1.csv", header = TRUE)

## Doing Coverage plots
# COI coverage
mat <- read.dna(file = "../alignments/COI653_65pc.fasta",
 format = "fasta", as.matrix = TRUE, as.character = TRUE)

for (i in 1:length(mat[1, ])){
	for (j in 1:length(mat[, 1])){
		if(mat[j, i] %in% c('a','t', 'g', 'c')){
			mat[j, i] <- "Base present"
		} else {
			mat[j, i] <- "Base absent"
		}
	}
}
cov_coi <- melt(mat)

# 65pc coverage
mat <- read.dna(file = "../alignments/65_percent_all_mitogenomes.gatk.nOG.aligned_CODES_sname 2.fasta", format = "fasta", as.matrix = T, as.character = T)

for (i in 1:length(mat[1, ])){
	for (j in 1:length(mat[, 1])){
		if(mat[j, i] %in% c('a','t', 'g', 'c')){
			mat[j, i] <- "Base present"
		} else {
			mat[j, i] <- "Base absent"
		}
	}
}
cov_65 <- melt(mat)


# 5pc coverage
#mat <- read.dna(file = "../alignments/5_percent_all_mitogenomes.gatk.nOG.aligned_sname 2.fasta", format = "fasta", as.matrix = T, as.character = T)
#mat <- read.dna(file = "../alignments/5_percent_all_mitogenomes.gatk.nOG.aligned_F_CODES_forR 2.fasta", format = "fasta", as.matrix = T, as.character = T)

for (i in 1:length(mat[1, ])){
	for (j in 1:length(mat[, 1])){
		if(mat[j, i] %in% c('a','t', 'g', 'c')){
			mat[j, i] <- "Base present"
		} else {
			mat[j, i] <- "Base absent"
		}
	}
}
cov_5 <- melt(mat)

# Panel figure with all of these

decay <- ggplot(dat) + geom_point(aes(x=Year.of.collection, y=Percent.Mitogenome.Coverage),
	pch=21, bg=alpha('steelblue', 0.5), size=3) + xlab('Year of collection') + ylab('Mitogenome coverage (%)') +
	ylim(0,100) + 
	theme_classic(base_size = 13) 

pc5 <- ggplot(cov_5, aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) + 
	xlab("Mitogenome Postion") + ylab("Individual") + 
	scale_fill_manual(values=c("#f0f0f0", "#636363"),name="") + scale_y_discrete(labels = NULL, breaks=NULL) +
	theme(panel.background = element_rect(fill = "white"), legend.position = "none")

pc65 <- ggplot(cov_65, aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) + 
	xlab("Mitogenome Postion") + ylab("Individual") + 
	scale_fill_manual(values=c("#f0f0f0", "#636363"),name="") + scale_y_discrete(labels = NULL, breaks=NULL) +
	theme(panel.background = element_rect(fill = "white"), legend.position = "none")

legend <- get_legend(ggplot(cov_coi, aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) + 
		xlab("COI Postion") + ylab("") + 
		scale_fill_manual(values=c("#f0f0f0", "#636363"),name="") + scale_y_discrete(labels = NULL, breaks=NULL) +
		theme(panel.background = element_rect(fill = "white")))

coi <- ggplot(cov_coi, aes(x = Var2, y = Var1)) + geom_raster(aes(fill=value)) + 
		xlab("COI Postion") + ylab("") + 
		scale_fill_manual(values=c("#f0f0f0", "#636363"),name="") + scale_y_discrete(labels = NULL, breaks=NULL) +
		theme(panel.background = element_rect(fill = "white"), legend.position = "none")


panel <- plot_grid(decay, pc5, pc65, coi, labels = "AUTO", nrow=2, ncol=2)

# Panel figure with all of these
tiff(file = "fig1.tiff", compression = "lzw")
	plot_grid(panel, legend, rel_widths = c(2, .4))
dev.off()
