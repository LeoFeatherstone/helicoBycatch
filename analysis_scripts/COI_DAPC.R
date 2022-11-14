
library(adegenet)

# PCA through Adegenet
COI653_Data = fasta2DNAbin("COI653_65pc.fasta")

GI_obj_COI_DAPC = DNAbin2genind(COI653_Data, polyThres = 0.01)
GI_obj_COI_DAPC

# Finding clusters in the data and denoting this as 'grp'. 
# We have 73 eigenvalues, so I will elect to keep 73 eigenvalues.
# It's not too taxing computationally, according to the tute since as we seek fewer than 10 clusters
# Each time choose 73 PC's and 2-7 clusters. Repeat this in the dapc funnction!

grp_COI_DAPC <- find.clusters(GI_obj_COI_DAPC, max.n.clust = 100, n.pca = 100)
# Seeing which cluster of (of choices 2-6) that each individual falls into
#write.csv(grp_COI_DAPC$grp, "COI_CLUSTER_ALLOCATION_SIX.csv")

# Doing the DAPC I'll then produce a few plots. I'll repeat this 
# for cluster numbers of 2 - 7

# Since 2 - 7 is a small number of clusters, it's ok to keep
# all the PCs (73)
dapc_COI <- dapc(GI_obj_COI_DAPC, grp_COI_DAPC$grp, n.pca = 73)
post <- dapc_COI$posterior
colnames(post) <- c("cluster1", "cluster2")

write.csv(dapc_COI$posterior, row.names = T, file = "COI_653_posteriorprob_k2.csv")
locs <- read.csv("~/Desktop/H.armigera project/MULTIVARIATE_ANALYSES/COI656_probs_locs.csv")


locsnew <- data.frame()

for (i in 1:length(locs$Sample)){
  
  if(any(grepl(locs$Sample[i], x = rownames(post)) > 0)){
    
    info <- c(i, locs$Sample[i], post[locs$Sample[i],], locs$Latitude[i], locs$Longitude[i])
    locsnew <- rbind(locsnew, info)
                     
  }
}
colnames(locsnew) <- c("i", "Sample", "cluster1", "cluster2", "Latitude", "Longitude")

locsnew <- cbind(locs[locsnew$i,]$Sample, locsnew)
locsnew <- locsnew[,c(1,4:7)]
colnames(locsnew) <- c("Sample", "cluster1", "cluster2", "Latitude", "Longitude")

write.csv(locsnew, file = "COI653_locs_TEMP.csv")
locsnew2 <- read.csv("COI653_locs_TEMP.csv")
