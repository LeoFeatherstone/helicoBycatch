############################################################################
############################# DACPC Figures ################################
############################################################################

library("oz")
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library("grid")

############### 5% DAPC ################
five_pc_aln <- fasta2DNAbin("./data/5percent_mitogenomes.fasta")
five_pc_genind <- DNAbin2genind(five_pc_aln, polyThres = 0.01)
five_pc_scalegen <- scaleGen(five_pc_genind, NA.method = "mean")
five_pc_pca <- dudi.pca(
  five_pc_scalegen,
  cent = FALSE,
  scale = FALSE,
  scannf = FALSE,
  nf = 3
) # We have 189 eigenvalues ==> 189 PC's

# Grouping. Set n.pca to no. eigenvalues (=189). Max Clusters = 10 to cover 6
five_pc_clustered <- find.clusters(
  five_pc_genind,
  max.n.clust = 10,
  n.pca = 189
) # choose two clusters

five_pc_dapc <- dapc(five_pc_genind, five_pc_clustered$grp, n.pca = 189)

mycol <- c("#4d4d4d", "#b2182b") # armigera red; conferta - grey

# test discriminant function
scatter(five_pc_dapc,
  scree.da = FALSE,
  bg = "white",
  pch = 20,
  cell = 0,
  cstar = 0,
  col = mycol,
  cex = 3,
  clab = 0,
  leg = TRUE,
  txt.leg = paste(c("Conferta", "Armigera"))
) # presents as a single discriminant function because k=2

# scatterplot
five_pc_scatterplot_data <- cbind.data.frame(
  five_pc_dapc$tab[, c(1, 2)],
  five_pc_dapc$grp
)
colnames(five_pc_scatterplot_data) <- c("PC1", "PC2", "Cluster")

five_pc_scatterplot <- ggplot(five_pc_scatterplot_data) +
  geom_point(aes(x = PC1, y = PC2, fill = Cluster), shape = 21, size = 5) +
  scale_fill_manual(values = alpha(mycol, 0.6), name = "") +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white")
  )


# 5pc map
five_pc_prob_locations <- read.csv("./data/5pc_mito_locations_probs.csv")

post <- five_pc_prob_locations[, 2:3]
post <- post[order(post[, 1]), ]

col <- vector()
for (i in seq_along(five_pc_prob_locations[, 1])) {
  if (five_pc_prob_locations$Subspecies[i] == "H. a. armigera") {
    col[i] <- "#b2182b"
  } else {
    col[i] <- "#4d4d4d"
  }
}

post <- five_pc_prob_locations[1:221, 2:3]
post <- post[, c(2, 1)]

layout_matrix <- matrix(c(1, 4, 2, 3), nrow = 2, ncol = 2)
layout(
  mat = layout_matrix,
  heights = c(1, 2),
  widths = c(2, 2)
)
layout.show(4)

plot.new()
legend(
  "center",
  legend = c(
    "Conferta",
    "Armigera"
  ),
  col = alpha(unique(col), 0.55),
  pch = 16,
  cex = 1.5,
  pt.cex = 2.5,
  border = FALSE
)


barplot(t(post[order(post[, 2]), ]),
  col = alpha(col, 0.55),
  border = NA,
  space = 0,
  xlab = "Individuals",
  ylab = "Posterior Probability",
  names.arg = rep(" ", length(post[, 1])),
  cex.lab = 1.25,
  mgp = c(1.5, 0, -0.9),
  las = 1
)

# Making Map
oz <- oz(states = TRUE)
points(
  jitter(five_pc_prob_locations$Lon_final, amount = 0.6),
  jitter(five_pc_prob_locations$Lat_final, amount = 0.6),
  cex = 2.5,
  pch = 16,
  col = alpha(col, 0.55)
)

vp_bottom_left <- viewport(
  height = unit(.5, "npc"),
  width = unit(0.5, "npc"),
  just = c("right", "top"),
  y = 0.5,
  x = 0.5
)

print(p5, vp = vp_bottom_left)


############### COI DAPC ################
coi_aln <- fasta2DNAbin("./data/COI_653bp_f.fasta")

coi_genind <- DNAbin2genind(coi_aln, polyThres = 0.01)
# Choose 73 PC's and 2-7 clusters. Repeat this in the dapc funnction

coi_clustered <- find.clusters(
  coi_genind,
  max.n.clust = 100,
  n.pca = 100
) # choose 2

coi_dapc <- dapc(coi_genind, coi_clustered$grp, n.pca = 73) # choose 10 DFs
post <- coi_dapc$posterior
colnames(post) <- c("cluster1", "cluster2") # cluster 1 = conferta; 2 = armigera

mycol <- c("#4d4d4d", "#b2182b") # armigera red; conferta - grey

# scatterplot
coi_scatterplot_data <- cbind.data.frame(coi_dapc$tab[, c(1, 2)], coi_dapc$grp)
colnames(coi_scatterplot_data) <- c("PC1", "PC2", "Cluster")

coi_scatterplot <- ggplot(coi_scatterplot_data) +
  geom_point(aes(x = PC1, y = PC2, fill = Cluster), shape = 21, size = 5) +
  scale_fill_manual(values = alpha(mycol, 0.6), name = "") +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white")
  )

coi_probs_locations <- read.csv("./data/COI653_locs_TEMP.csv")
post <- coi_probs_locations[, 3:4]
post <- post[order(post[, 1]), ]

col <- vector()
for (i in seq_along(coi_probs_locations$Subspecies)) {
  if (coi_probs_locations$Subspecies[i] == "H. a. armigera") {
    col[i] <- "#b2182b"
  } else {
    col[i] <- "#4d4d4d"
  }
}

layout_matrix <- matrix(c(1, 4, 2, 3), nrow = 2, ncol = 2)
layout(
  mat = layout_matrix,
  heights = c(1, 2),
  widths = c(2, 2)
)
layout.show(4)

plot.new()
legend(
  "center",
  legend = c(
    "Conferta",
    "Armigera"
  ),
  col = alpha(unique(col), 0.55),
  pch = 16,
  cex = 1.5,
  pt.cex = 2.5,
  border = FALSE
)


barplot(
  t(post[order(post[, 2]), ]),
  col = alpha(col, 0.55),
  border = NA,
  space = 0,
  xlab = "Individuals",
  ylab = "Posterior Probability",
  names.arg = rep(" ", length(post[, 1])),
  cex.lab = 1.25,
  mgp = c(1.5, 0, -0.9),
  las = 1
)


# Making Map
oz <- oz(states = TRUE)
oz <- points(
  jitter(coi_probs_locations$Longitude, amount = 0.6),
  jitter(coi_probs_locations$Latitude, amount = 0.6),
  cex = 2.5,
  pch = 16,
  col = alpha(col, 0.55)
)

vp_bottom_left <- viewport(
  height = unit(.5, "npc"),
  width = unit(0.5, "npc"),
  just = c("right", "top"),
  y = 0.5,
  x = 0.5
)

print(coi_scatterplot, vp = vp_bottom_left)

dev.off()


################### 65% DAPC #########################

sixtyfive_pc_aln <- fasta2DNAbin("./data/65percent_mitogenomes.fasta")
sixtyfive_pc_genind <- DNAbin2genind(sixtyfive_pc_aln, polyThres = 0.01)


sixtyfive_pc_scalegen <- scaleGen(sixtyfive_pc_genind, NA.method = "mean")

sixtyfive_pc_pca <- dudi.pca(
  sixtyfive_pc_scalegen,
  cent = FALSE,
  scale = FALSE,
  scannf = FALSE,
  nf = 3
)
# We have 54 eigenvalues ==> 54 PC's
# Set num pca to num eigenvalues (=54). Max Clusters = 10 to cover 6
sixtyfive_pc_clustered <- find.clusters(
  sixtyfive_pc_genind,
  max.n.clust = 10,
  n.pca = 54
) # choose two clusters

# Doing the DAPC
sixtyfive_pc_dapc <- dapc(
  sixtyfive_pc_genind,
  sixtyfive_pc_clustered$grp,
  n.pca = 54
) # choose 10 DFs

mycol2 <- c("#4d4d4d", "#b2182b") # armigera red; conferta - grey

# scatter
sixtyfive_pc_scatterplot_data <- cbind.data.frame(
  sixtyfive_pc_dapc$tab[, c(1, 2)],
  sixtyfive_pc_dapc$grp
)
colnames(sixtyfive_pc_scatterplot_data) <- c("PC1", "PC2", "Cluster")

sixtyfive_pc_scatterplot <- ggplot(sixtyfive_pc_scatterplot_data) +
  geom_point(aes(x = PC1, y = PC2, fill = Cluster), shape = 21, size = 5) +
  scale_fill_manual(values = alpha(mycol, 0.6), name = "") +
  theme(
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white")
  )

# Mapping 65% Data for preprint
sixtyfive_probs_locs <- read.csv("./data/65pc_dapc_locations.csv")

post <- sixtyfive_probs_locs[, c(4, 5)]
post <- post[order(post[, 1]), ]

col <- vector()
for (i in seq_along(sixtyfive_probs_locs$Subspecies)) {
  if (sixtyfive_probs_locs$Subspecies[i] == "H. a. armigera") {
    col[i] <- "#b2182b"
  } else {
    col[i] <- "#4d4d4d"
  }
}

layout_matrix <- matrix(c(1, 4, 2, 3), nrow = 2, ncol = 2)
layout(
  mat = layout_matrix,
  heights = c(1, 2),
  widths = c(2, 2)
)
layout.show(4)

plot.new()
legend(
  "center",
  legend = c(
    "Conferta",
    "Armigera"
  ),
  col = alpha(unique(col), 0.55),
  pch = 16,
  cex = 1.5,
  pt.cex = 2.5,
  border = FALSE
)


barplot(
  t(post[order(post[, 2]), ]),
  col = alpha(col, 0.55),
  border = NA,
  space = 0,
  xlab = "Individuals",
  ylab = "Posterior Probability",
  names.arg = rep(" ", length(post[, 1])),
  cex.lab = 1.25,
  mgp = c(1.5, 0, -0.9),
  las = 1
)

# Making Map
oz <- oz(states = TRUE)
points(jitter(sixtyfive_probs_locs$Longitude, amount = 0.6),
  jitter(sixtyfive_probs_locs$Latitude, amount = 0.6),
  cex = 2.5,
  pch = 16,
  col = alpha(col, 0.55)
)

vp_bottom_left <- viewport(
  height = unit(.5, "npc"),
  width = unit(0.5, "npc"),
  just = c("right", "top"),
  y = 0.5,
  x = 0.5
)

print(sixtyfive_pc_scatterplot, vp = vp_bottom_left)

dev.off()
