### A script to do the DAPC for the >/= 65% coverage data. No-read Bias analysis first ###

#install.packages("adegenet", dep=TRUE)
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")


setwd("~/Desktop/H.armigera project/MULTIVARIATE_ANALYSES/")

################################################
#### Doing the bias check for missing data ####
#sixty_five_pc_dna_mat <- read.dna(format = "fasta", as.matrix = T, as.character = T,
                                  file = "65_percent_all_mitogenomes.gatk.nOG.aligned_CODES.fasta")
#class(sixty_five_pc_dna_mat)
#write.table(sixty_five_pc_dna_mat, file = "sixtyfive_pc_mito_mat_temp.txt")

binary_sixtyfive_pc_mat <- as.matrix(read.table("sixtyfive_pc_mito_mat_temp.txt"))
class(binary_sixtyfive_pc_mat) <- "numeric"


dim(binary_sixtyfive_pc_mat)

# Doing the PCA
sixtyfive_pc_pca <- prcomp(binary_sixtyfive_pc_mat)

# Calculating the proportion of covered sites for each individual. Store in vector.
sixtyfive_pc_prop_covered <- vector()

for (i in 1:length(binary_sixtyfive_pc_mat[,1])){
  sixtyfive_pc_prop_covered[i] <- sum(binary_sixtyfive_pc_mat[i,])/length(binary_sixtyfive_pc_mat[1,])
}
length(sixtyfive_pc_prop_covered)

# Now making a data frame to present this with a colout scale 

sixtyfive_pc_QC_df <- as.data.frame(cbind(sixtyfive_pc_pca$x[,1], sixtyfive_pc_pca$x[,2], sixtyfive_pc_pca$x[,3], sixtyfive_pc_pca$x[,4], sixtyfive_pc_pca$x[,5],
                                          sixtyfive_pc_pca$x[,6], sixtyfive_pc_prop_covered))
var_exp_65 <- 100*sixtyfive_pc_pca$sdev^2/sum(sixtyfive_pc_pca$sdev^2)
var_exp_65


names(sixtyfive_pc_QC_df) <- c("PC1", "PC2","PC3", "PC4", "PC5", "PC6", "Coverage")
#names(COI_QC_df) <- c("PC1 59.2%", "PC2 18.5%","PC3 5.11%", "PC4 2.80%", "PC5 1.36%", "PC6 1.15%", 
"PC7 0.098%", "PC8 0.075%", "Coverage")

sixtyfive_pc_coverage_plot_1 <- ggplot(sixtyfive_pc_QC_df, aes(x = PC1, y = PC2, colour = Coverage)) + geom_point()+
  xlab("PC1 - 37%") + ylab("PC2 - 5%") + guides(fill = guide_colourbar(title = "Coverage", ticks = F, title.position = "top")) +
  theme(axis.title.x =element_text(size=15), axis.title.y =element_text(size=15), legend.title = element_text(size=15),
        axis.text.y = element_text(size=13), axis.text.x = element_text(size=13))
sixtyfive_pc_coverage_plot_1
sixtyfive_pc_coverage_plot_2 <- ggplot(sixtyfive_pc_QC_df, aes(x = PC3, y = PC2, colour = Coverage)) + geom_point()+
  xlab("PC2 - 5%") + ylab("PC3 - 4%") + guides(fill = guide_colourbar(title = "Coverage", ticks = F, title.position = "top")) +
  theme(axis.title.x =element_text(size=15), axis.title.y =element_text(size=15), legend.title = element_text(size=15),
        axis.text.y = element_text(size=13), axis.text.x = element_text(size=13))
sixtyfive_pc_coverage_plot_2
sixtyfive_pc_coverage_plot_3 <- ggplot(sixtyfive_pc_QC_df, aes(x = PC4, y = PC3, colour = sixtyfive_pc_prop_covered)) + geom_point()
sixtyfive_pc_coverage_plot_3
sixtyfive_pc_coverage_plot_4 <- ggplot(sixtyfive_pc_QC_df, aes(x = PC5, y = PC6, colour = sixtyfive_pc_prop_covered)) + geom_point()
sixtyfive_pc_coverage_plot_4

#install.packages("gridExtra")
#library(gridExtra)

grid.arrange(sixtyfive_pc_coverage_plot_1, sixtyfive_pc_coverage_plot_2)

################################
### Now looking at gene coverage 

par(mfrow = c(1,1))
hist(sixtyfive_pc_prop_covered, freq = T)

# Now making a plot of coverage across the gene 
sixtyfive_pc_num_ind_site <- vector()
length(binary_sixtyfive_pc_mat[1,])
for (i in 1:length(binary_sixtyfive_pc_mat[1,])){
  sixtyfive_pc_num_ind_site[i] <- sum(binary_sixtyfive_pc_mat[,i])
}
length(sixtyfive_pc_num_ind_site)
sites <- seq(1,length(sixtyfive_pc_num_ind_site), by = 1)
length(num_ind_site)
sixtyfive_pc_site_cov_df <- as.matrix(cbind(sixtyfive_pc_num_ind_site, as.factor(sites)))

par(mfrow = c(1,1))
barplot(sixtyfive_pc_site_cov_df[,1], ylab = "Number of ndividuals covering mitogenome position", xlab = "Mitogenome")
# Subsetting desirable sites
sixtyfive_pc_enough_reads <- sixtyfive_pc_site_cov_df[sixtyfive_pc_site_cov_df[,1] > 350,]
sixtyfive_pc_enough_reads[501:540,]
sixtyview(five_pc_enough_reads)
### Sites to keep therefore are ###

############################################################################################################
### Doing the DAPC ###


SIXTYFIVE_PERCENT_MITO_Data = fasta2DNAbin("65_percent_all_mitogenomes.gatk.nOG.aligned_CODES.fasta")
SIXTYFIVE_PERCENT_MITO_GIO = DNAbin2genind(SIXTYFIVE_PERCENT_MITO_Data, polyThres = 0.01)


SIXTYFIVE_PERCENT_MITO_PCA_OBJ = scaleGen(SIXTYFIVE_PERCENT_MITO_GIO, NA.method = "mean")

SIXTYFIVE_PERCENT_MITO_PCA <- dudi.pca(SIXTYFIVE_PERCENT_MITO_PCA_OBJ, cent = FALSE, scale = FALSE, scannf = FALSE, nf=3)
SIXTYFIVE_PERCENT_MITO_PCA

#########################################################################################################
############################################# Doing my own PCA Plots ####################################
#########################################################################################################


pca_col_65 <- vector()

for (i in 1:length(dapc_SIXTYFIVE_PERCENT_MITO$assign)){
  if (dapc_SIXTYFIVE_PERCENT_MITO$assign[i] == 1){
  pca_col_65[i] <- as.character(expression(italic("H. a. armigera")))
  } else {
    pca_col_65[i] <- as.character(expression(italic("H. a. conferta")))
  }
}

# Setting colour Data
Sixtyfive_PCA_df <- cbind(SIXTYFIVE_PERCENT_MITO_PCA$li, dapc_SIXTYFIVE_PERCENT_MITO$assign)
write.csv(Sixtyfive_PCA_df, file = "65PC_SelfPCA_PLOT.csv")
Sixtyfive_PCA_df <- read.csv("65PC_SelfPCA_PLOT.csv")

par(mgp = c(2.75,1,0))
plot(x=Sixtyfive_PCA_df$Axis1, y=Sixtyfive_PCA_df$Axis2, col = "black", bg = as.vector(Sixtyfive_PCA_df$SubSpCOL),
     xlab = "PC1 - 8%", ylab = "PC2 = 6%", cex = 1.5, pch = 21, cex.axis = 1.2, cex.lab = 1.5, las = 1)

legend("bottomright", legend = c(expression(italic("H. a. armigera")), expression(italic("H. a. conferta"))),
       pch = 21, col = "black", pt.bg = c("#2ca25f", "#43a2ca"), cex = 1.5)

plot(x=Sixtyfive_PCA_df$Axis2, y=Sixtyfive_PCA_df$Axis3, col = "black", bg = as.vector(Sixtyfive_PCA_df$SubSpCOL),
     xlab = "PC2 - 6%", ylab = "PC3 = 5%", cex = 1.5, pch = 21, cex.axis = 1.2, cex.lab = 1.5, las = 1)

legend("bottomright", legend = c(expression(italic("H. a. armigera")), expression(italic("H. a. conferta"))),
       pch = 21, col = "black", pt.bg = c("#2ca25f", "#43a2ca"), cex = 1.5)

SIXTYFIVE_PERCENT_MITO_PCA

##### Getting outliers: ######
# 65%: under 10 on PC1 and above -10 on PC2
out_65 <- Sixtyfive_PCA_df[!((Sixtyfive_PCA_df$Axis1 < 10) & (Sixtyfive_PCA_df$Axis2 > -10)),]
write.csv(out_65, file = "65pc_PCA_OUTLIERS.csv")

# We have 139 eigenvalues ==> 139 PC's

# Now to do the DAPC (RETURN LATER TO GET FIGURES FROM THE PCA)

# Grouping. Set n.pca to no. eigenvalues (=139). Max Clusters = 10 to cover 6

grp_SIXTYFIVE_PERCENT_MITO <- find.clusters(SIXTYFIVE_PERCENT_MITO_GIO, max.n.clust = 10, n.pca = 139)
# Doing the DAPC
dapc_SIXTYFIVE_PERCENT_MITO <- dapc(SIXTYFIVE_PERCENT_MITO_GIO, grp_SIXTYFIVE_PERCENT_MITO$grp, n.pca = 139)

dapc_SIXTYFIVE_PERCENT_MITO$posterior
length(binary_sixtyfive_pc_mat[1,])
dim(binary_sixtyfive_pc_mat)

order_65 <- cbind(binary_sixtyfive_pc_mat, dapc_SIXTYFIVE_PERCENT_MITO$assign)
dim(order_65)
order_65[,15438]
binary_sixtyfive_pc_mat_ord <- order_65[order(order_65[,15438]),]
dim(binary_sixtyfive_pc_mat_ord)
binary_sixtyfive_pc_mat_ord[,15438]
# Removing order column
binary_sixtyfive_pc_mat_ord <- binary_sixtyfive_pc_mat_ord[,-15438]
cols_65 <- c(t(binary_sixtyfive_pc_mat_ord))

indiv_65 <- vector()
site_65 <- vector()

for(i in 1:length(binary_sixtyfive_pc_mat_ord[,1])){
  site_65 <- c(site_65, seq(1:length(binary_sixtyfive_pc_mat_ord[1,])))
}

for(i in 1:length(binary_sixtyfive_pc_mat_ord[,1])){
  indiv_65 <- c(indiv_65, rep(i, length = length(binary_sixtyfive_pc_mat_ord[1,])))
}

length(cols_65)
point_data_65 <- as.data.frame(cbind(indiv_65, site_65, cols_65))

## Doing levelplot
library(lattice)
levelplot(cols_65 ~ site_65*indiv_65, data = point_data_65, colorkey = F, col.regions = c(5,6), labels = c("Read", "No Read"),
          xlab = list("Position", cex = 1.5), ylab = list("Individual Moths", cex = 1.5), 
          scales = list(x = list(cex = 1.2), y = list(cex = 1.2), tck = c(1,0)) )

levelplot(cols ~ site*indiv, data = point_data, colorkey = F,
          col.regions = c(6,5), labels = c("Read", "No Read"),
          xlab = list("Position in CO1", cex = 1.5), ylab = list("Individual Moths", cex = 1.5), 
          scales = list(x = list(cex = 1.2), y = list(cex = 1.2), tck = c(1,0)),)

abline(h = 31)
length(dapc_SIXTYFIVE_PERCENT_MITO[dapc_SIXTYFIVE_PERCENT_MITO$assign == 1])

write.csv(cbind(dapc_SIXTYFIVE_PERCENT_MITO$posterior, dapc_SIXTYFIVE_PERCENT_MITO$assign), file = "65pc_dapc_locations.csv")

write.csv(cbind(SIXTYFIVE_PERCENT_MITO_PCA$li$Axis1), file = "65_PC1_values.csv")

################################################################################################################################################
################################################ Doing the DAPC Plot ###########################################################################
################################################################################################################################################

# Run DAPC for 656 first and then do barplot 
DAPC_65_post_mat <- as.matrix(dapc_SIXTYFIVE_PERCENT_MITO$posterior)
DAPC_65_post_mat <- DAPC_65_post_mat[order(DAPC_65_post_mat[,1]),]

par(oma = c(0,1,0,0), mgp = c(3,1,0))
barplot(t(DAPC_65_post_mat), col = c("#2ca25f", "#43a2ca"), border = NA, space = 0.0, 
        names.arg = rep("", length(t(DAPC_65_post_mat)[1,])), cex.axis = 1.2, las = 1)
title(xlab = "Individual Moths", cex.lab = 1.5, line = 1)
title(ylab = "Posterior Probability of Cluster Membership", cex.lab = 1.5, line = 3)
par(oma = c(0,1,1,2))
legend("topright", legend = c(expression(italic("H. a. armigera")), expression(italic("H. a. conferta"))), pch = 15,
       col = c("#43a2ca" , "#2ca25f"), cex = 1.2, bty = "o", box.col = "white")
par(oma = c(0,0,0,0), mgp = c(3,1,0))      
        