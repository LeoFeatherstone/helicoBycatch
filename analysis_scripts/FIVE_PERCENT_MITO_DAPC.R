##### 5% MITOGENOME PCA. CLUSTERS 2 - 6 ##########

# Doing a PCA of the 5% mitogenome data to find how many PC's there are

#install.packages("adegenet", dep=TRUE)
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")


setwd("~/Desktop/H.armigera project/MULTIVARIATE_ANALYSES/")

FIVE_PERCENT_MITO_Data = fasta2DNAbin("5_percent_all_mitogenomes.gatk.nOG.aligned_F_CODES_forR.fasta")
FIVE_PERCENT_MITO_GIO = DNAbin2genind(FIVE_PERCENT_MITO_Data, polyThres = 0.01)


FIVE_PERCENT_MITO_PCA_OBJ = scaleGen(FIVE_PERCENT_MITO_GIO, NA.method = "mean")

FIVE_PERCENT_MITO_PCA <- dudi.pca(FIVE_PERCENT_MITO_PCA_OBJ, cent = FALSE, scale = FALSE, scannf = FALSE, nf=3)
FIVE_PERCENT_MITO_PCA
# We have 139 eigenvalues ==> 139 PC's

# Now to do the DAPC (RETURN LATER TO GET FIGURES FROM THE PCA)

# Grouping. Set n.pca to no. eigenvalues (=139). Max Clusters = 10 to cover 6

grp_FIVE_PERCENT_MITO <- find.clusters(FIVE_PERCENT_MITO_GIO, max.n.clust = 10, n.pca = 139)
write.csv(grp_FIVE_PERCENT_MITO$grp, "FIVE_PERCENT_MITO_TWO_CLUSTER_ALLOCATION.csv")
# Doing the DAPC
dapc_FIVE_PERCENT_MITO <- dapc(FIVE_PERCENT_MITO_GIO, grp_FIVE_PERCENT_MITO$grp, n.pca = 139)

mycol = funky(2)
scatter.dapc(dapc_FIVE_PERCENT_MITO, scree.da = FALSE, bg="white", pch = 20, cell = 0, cstar = 0,
             col = mycol, cex=3,clab=0,leg=TRUE,txt.leg = paste("Cluster", 1:2))

#############################################################################################################
########################### Doing my own PCA plots now ######################################################
#############################################################################################################
Five_PCA_df <- cbind(FIVE_PERCENT_MITO_PCA$li, dapc_FIVE_PERCENT_MITO$assign)
write.csv(Five_PCA_df, file = "5PC_SelfPCA_PLOT.csv")
Five_PCA_df <- read.csv("5PC_SelfPCA_PLOT.csv")

par(mgp = c(2.75,1,0))
plot(x=Five_PCA_df$Axis1, y=Five_PCA_df$Axis2, col = "black", bg = as.vector(Five_PCA_df$SubSpCol),
     xlab = "PC1 - 8%", ylab = "PC2 = 7%", cex = 1.5, pch = 21, cex.axis = 1.2, cex.lab = 1.5, las = 1)

legend("bottomright", legend = c(expression(italic("H. a. armigera")), expression(italic("H. a. conferta"))),
       pch = 21, col = "black", pt.bg = c("#2ca25f", "#43a2ca"), cex = 1.5)

plot(x=Five_PCA_df$Axis2, y=Five_PCA_df$Axis3, col = "black", bg = as.vector(Five_PCA_df$SubSpCol),
     xlab = "PC2 - 7%", ylab = "PC3 = 4%", cex = 1.5, pch = 21, cex.axis = 1.2, cex.lab = 1.5, las = 1)

legend("topleft", legend = c(expression(italic("H. a. armigera")), expression(italic("H. a. conferta"))),
       pch = 21, col = "black", pt.bg = c("#2ca25f", "#43a2ca"), cex = 1.5)

#### Gettong PCA outliers #####
# 5%: under 20 on PC1 and above -10 on PC2

out_5 <- Five_PCA_df[!((Five_PCA_df$Axis1 < 20) & (Five_PCA_df$Axis2 > -10)),]
out_5
write.csv(out_5, file = "5pc_PCA_OUTLIERS.csv")
################################################
#### Doing the bias check for missing data ####

five_pc_dna_mat <- read.dna(format = "fasta", as.matrix = T, as.character = T, file = "5_percent_all_mitogenomes.gatk.OG.aligned_F_CODES_forR.fasta")
class(five_pc_dna_mat)
write.table(five_pc_dna_mat, file = "five_pc_mito_mat_temp.txt")

binary_five_pc_mat <- as.matrix(read.table("five_pc_mito_mat_temp.txt"))
class(binary_five_pc_mat) <- "numeric"


dim(binary_five_pc_mat)

# Doing the PCA
five_pc_pca <- prcomp(binary_five_pc_mat)

# Calculating the proportion of covered sites for each individual. Store in vector.
five_pc_prop_covered <- vector()

for (i in 1:length(binary_five_pc_mat[,1])){
  five_pc_prop_covered[i] <- sum(binary_five_pc_mat[i,])/length(binary_five_pc_mat[1,])
}
length(five_pc_prop_covered)

# Now making a data frame to present this with a colout scale 

five_pc_QC_df <- as.data.frame(cbind(five_pc_pca$x[,1], five_pc_pca$x[,2],five_pc_pca$x[,3],five_pc_pca$x[,4],five_pc_pca$x[,5],
                                     five_pc_pca$x[,6], five_pc_prop_covered))
var_exp_5 <- 100*five_pc_pca$sdev^2/sum(five_pc_pca$sdev^2)
var_exp_5


names(five_pc_QC_df) <- c("PC1", "PC2","PC3", "PC4", "PC5", "PC6", "Coverage")
#names(COI_QC_df) <- c("PC1 59.2%", "PC2 18.5%","PC3 5.11%", "PC4 2.80%", "PC5 1.36%", "PC6 1.15%", 
"PC7 0.098%", "PC8 0.075%", "prop_covered")

five_pc_coverage_plot_1 <- ggplot(five_pc_QC_df, aes(x = PC1, y = PC2, colour = Coverage)) + geom_point()+
  xlab("PC1 - 53%") + ylab("PC2 - 4%") + guides(fill = guide_colourbar(title = "Coverage", ticks = F, title.position = "top")) +
  theme(axis.title.x =element_text(size=15), axis.title.y =element_text(size=15), legend.title = element_text(size=15),
        axis.text.y = element_text(size=13), axis.text.x = element_text(size=13))
five_pc_coverage_plot_1
five_pc_coverage_plot_2 <- ggplot(five_pc_QC_df, aes(x = PC2, y = PC3, colour = Coverage)) + geom_point()+
  xlab("PC2 - 4%") + ylab("PC3 - 3%") + guides(fill = guide_colourbar(title = "Coverage", ticks = F, title.position = "top")) +
  theme(axis.title.x =element_text(size=15), axis.title.y =element_text(size=15), legend.title = element_text(size=15),
        axis.text.y = element_text(size=13), axis.text.x = element_text(size=13))
five_pc_coverage_plot_2

five_pc_coverage_plot_3 <- ggplot(five_pc_QC_df, aes(x = PC4, y = PC3, colour = prop_covered)) + geom_point()
five_pc_coverage_plot_3
five_pc_coverage_plot_4 <- ggplot(five_pc_QC_df, aes(x = PC5, y = PC6, colour = prop_covered)) + geom_point()
five_pc_coverage_plot_4

#install.packages("gridExtra")
#library(gridExtra)

grid.arrange(five_pc_coverage_plot_1, five_pc_coverage_plot_2)

################################
### Now looking at gene coverage 

par(mfrow = c(1,1))
hist(five_pc_prop_covered, freq = T)

# Now making a plot of coverage across the gene 
five_pc_num_ind_site <- vector()
length(binary_five_pc_mat[1,])
for (i in 1:length(binary_five_pc_mat[1,])){
  five_pc_num_ind_site[i] <- sum(binary_five_pc_mat[,i])
}
length(five_pc_num_ind_site)
sites <- seq(1,length(five_pc_num_ind_site), by = 1)
length(num_ind_site)
five_pc_site_cov_df <- as.matrix(cbind(five_pc_num_ind_site, as.factor(sites)))

par(mfrow = c(1,1))
barplot(five_pc_site_cov_df[,1], ylab = "Number of ndividuals covering mitogenome position", xlab = "Mitogenome")
# Subsetting desirable sites
five_pc_enough_reads <- five_pc_site_cov_df[five_pc_site_cov_df[,1] > 350,]
five_pc_enough_reads[501:540,]
view(five_pc_enough_reads)
### Sites to keep therefore are ###

##### Standard 2 cluster DAPC plot first #####
mito5_probs <- (dapc_FIVE_PERCENT_MITO$posterior)
## Above object is a matrix 
mito5_probs_ordered <- mito5_probs[order(mito5_probs[,1], decreasing = T),]
mito5_probs_ordered %>% head()
mito5_probs_ordered_T <- t(mito5_probs_ordered)
mito5_probs_ordered_T
length(mito5_probs_ordered_T[1,])
par(oma = c(3.5,0,0,0))
barplot(mito5_probs_ordered_T, beside = F, ylab = "Posterior Probability of Cluster Assignment", xlab = "Individual Moths", 
        col = c(4,3), width = 1, names.arg = rep("", 256), space = 0,  border = NA, legend = F)
par(oma = c(0,0,0,0))
legend("bottomright" ,legend = c("H. armigera", "H. Conferta"), pch = 15, col = c(4,3), horiz = F)


##### Extracting admixed individuals 
mito5_probs_ordered_admixed <- mito5_probs_ordered[((mito5_probs_ordered[,1] > 0.01) & (mito5_probs_ordered[,1] < 0.99)),]
mito5_probs_ordered_admixed
write.csv(mito5_probs_ordered_admixed, file = "mito5_admixed.csv")

# Now getting one of those membership probability plots
compoplot(dapc_FIVE_PERCENT_MITO, txt.leg=paste("Cluster", 1:2), 
          lab="", xlab="individuals", col = mycol)

# Seeing probability distributions
hist(dapc_FIVE_PERCENT_MITO$posterior)


##### Doing a coverage matrix #######

# Start with binary matrix:
dim(binary_five_pc_mat)

indiv <- vector()
site <- vector()
cols <- vector()

dim(binary_five_pc_mat)
cols <- c(t(binary_five_pc_mat))

for (i in 1:256){
  indiv <- c(indiv, rep(i, length = 15438))
}

for (j in 1:256){
  len <- c(1:15438)
  site <- c(site, len)
}

length(cols)

# Trying a levelplot
# Ordering the matrix first so that I get groupings:
dapc_FIVE_PERCENT_MITO$assign
binary_five_pc_mat[,1]

ordering <- cbind(binary_five_pc_mat, dapc_FIVE_PERCENT_MITO$assign)
dim(ordering)
binary_five_pc_mat_ord <- ordering[order(ordering[,15439]),]
dim(binary_five_pc_mat_ord)
binary_five_pc_mat_ord[,15439]

# now removing order columns and re-assigning ordered individuals
binary_five_pc_mat_ord <- binary_five_pc_mat_ord[,-15439]
cols <- c(t(binary_five_pc_mat_ord))

library(lattice)

point_data <- as.data.frame(cbind(indiv, site, cols))

levelplot(cols ~ site*indiv, data = point_data, colorkey = F,
          col.regions = c(5,6),
          xlab = list("Position ", cex = 1.5), ylab = list("Individual Moths", cex = 1.5), 
          scales = list(x = list(cex = 1.2), y = list(cex = 1.2), tck = c(1,0)),)


abline(h = 201)


length(dapc_FIVE_PERCENT_MITO[dapc_FIVE_PERCENT_MITO$assign == 1])

###################################################################################
###################################################################################
####################### raster vis coverage plot ##################################
###################################################################################
###################################################################################

r <- raster(x = binary_five_pc_mat_ord, xmn=0, xmx=length(binary_five_pc_mat_ord[1,]), 
ymn=0, ymx=length(binary_five_pc_mat_ord[,1]), crs=NA, template=NULL)

r <- ratify(r)
rat <- levels(r)[[1]]

rat$landcover <- c("No Read", "Read")
levels(r) <- rat
extent <- c(0, length(binary_five_pc_mat_ord[1,]), 0, length(binary_five_pc_mat_ord[,1]))

dev.off()
levelplot(r, col.regions = c(5,6), axis.margin = T,
          xlab = list("Position in CO1", cex = 1.5), 
          ylab = list("Individual Moths", cex = 1.5),
          scales = list(x = list(cex = 1.2, lim = c(1:length(binary_five_pc_mat_ord[1,]))), 
                        y = list(cex = 1.2, lim = c(1:length(binary_five_pc_mat_ord[,1]))),
                        colorkey=list(labels=list(cex=1.2)))) 


################################################################################################################################################
################################################ Doing the DAPC Plot ###########################################################################
################################################################################################################################################

# Run DAPC for 656 first and then do barplot 
DAPC_5_post_mat <- as.matrix(dapc_FIVE_PERCENT_MITO$posterior)
DAPC_5_post_mat <- DAPC_5_post_mat[order(DAPC_5_post_mat[,2]),]

par(oma = c(0,1,0,0), mgp = c(3,1,0))
barplot(t(DAPC_5_post_mat), col = c("#43a2ca","#2ca25f"), border = NA, space = 0.0, 
        names.arg = rep("", length(t(DAPC_5_post_mat)[1,])), cex.axis = 1.2, las = 1)
title(xlab = "Individual Moths", cex.lab = 1.5, line = 1)
title(ylab = "Posterior Probability of Cluster Membership", cex.lab = 1.5, line = 3)
par(oma = c(0,1,1,2))
legend("topright", legend = c(expression(italic("H. a. armigera")), expression(italic("H. a. conferta"))), pch = 15,
       col = c("#43a2ca" , "#2ca25f"), cex = 1.2, bty = "o", box.col = "white")
par(oma = c(0,0,0,0), mgp = c(3,1,0))      



