## Plotting for figures 2 and sups
library(ggtree)
library(ggplot2)
library(treeio)

setwd(paste0(getwd(), "/figures"))

# 5pc
five_pc <- treeio::read.beast("5percent_mitogenomes_relaxed_MCC.trees")

# removing negative branch lengths??

# get locations for tips
loc <- vector()
for (i in 1:length(five_pc@phylo$tip.label)){
  if (grepl("AUS_", five_pc@phylo$tip.label[i]) || grepl("NZL_", five_pc@phylo$tip.label[i])){
    loc[i] <- "Australian"
  } else {
    loc[i] <- "Rest of World"
  }
}

loc <- as.data.frame(loc)
colnames(loc) <- c("Location")
rownames(loc) <- five_pc@phylo$tip.label

# get colour vector for tips
col <- vector()
for (i in 1:length(loc$Location)){
  if(grepl("AUS", loc$Location[i]) || grepl("NZL", loc$Location[i])){
    col[i] <- "#b2182b"
  } else {
    col[i] <- "#4d4d4d"
  }
}

# setting up basic tree plot
p5 <- ggtree(five_pc, aes(color = posterior)) +
 scale_color_continuous(name = "Posterior\nProbability") +
 #geom_cladelabel(node = 407, label = "Rest of \n World", offset = 15, offset.text = 15, fontsize = 5,
  #col = "#b2182b") +
 #geom_cladelabel(node = 261, label = "Australian", offset = 15, offset.text = 15, fontsize = 5,
  #col  = "#4d4d4d") +
 geom_treescale(x = 0, y = 100, label = "yrs") +
 theme(legend.position = "left")


tiff(file = "figS2.tiff", compression = "lzw")
gheatmap(p5, data = loc, width = 0.2, offset = 650,
         colnames = F) +
  scale_fill_manual(values = c("#4d4d4d",  "#b2182b"), name = "") +
  theme(legend.text = element_text(size = 12)) +
  scale_x_ggtree(breaks = c(0, 1000, 2000)) 
dev.off()

# 65 pc
sixtyfive_pc <- treeio::read.beast("65percent_mitogenomes_relaxed_MCC.trees")
# removing negative branch lengths??

# get locations for tips
loc <- vector()
for (i in 1:length(sixtyfive_pc@phylo$tip.label)){
  if (grepl("AUS_", sixtyfive_pc@phylo$tip.label[i]) || grepl("NZL_", sixtyfive_pc@phylo$tip.label[i])){
    loc[i] <- "Australian"
  } else {
    loc[i] <- "Rest of World"
  }
}

loc <- as.data.frame(loc)
colnames(loc) <- c("Location")
rownames(loc) <- sixtyfive_pc@phylo$tip.label

# get colour vector for tips
col <- vector()
for (i in 1:length(loc$Location)){
  if(grepl("AUS", loc$Location[i]) || grepl("NZL", loc$Location[i])){
    col[i] <- "#b2182b"
  } else {
    col[i] <- "#4d4d4d"
  }
}

# setting up basic tree plot
p65 <- ggtree(sixtyfive_pc, aes(color = posterior)) +
 scale_color_continuous(name = "Posterior\nProbability") +
 #geom_cladelabel(node = 407, label = "Rest of \n World", offset = 15, offset.text = 15, fontsize = 5,
  #col = "#b2182b") +
 #geom_cladelabel(node = 261, label = "Australian", offset = 15, offset.text = 15, fontsize = 5,
  #col  = "#4d4d4d") +
 geom_treescale(x = 0, y = 50, label = "yrs") +
 theme(legend.position = "left") 
 #geom_text(aes(label=node), hjust=-.3) +
  #geom_cladelabel(node = 57, label = "Rest of \n World", offset = 15, offset.text = 15, fontsize = 5,
        #          col = "#b2182b") #+
  #geom_cladelabel(node = 71, label = "Australian", offset = 15, offset.text = 15, fontsize = 5,
  # col  = "#4d4d4d")


tiff(file = "fig3.tiff", compression = "lzw")
gheatmap(p65, data = loc, width = 0.2, offset = 50,
         colnames = FALSE) +
  scale_fill_manual(values = c("#4d4d4d",  "#b2182b"), name = "") +
  theme(legend.text = element_text(size = 12),
  legend.title = element_text(size = 12)) +
  scale_x_ggtree(breaks = seq(from = 1600, to = 2000, by = 100)) 
dev.off()


sixtyfive_pc <- read.nexus("65%_relaxed_HPD.tree")

# get locations for tips
loc <- vector()
for (i in 1:length(sixtyfive_pc$tip.label)){
  if (grepl("AUS_", sixtyfive_pc$tip.label[i]) || grepl("NZL_", sixtyfive_pc$tip.label[i])){
    loc[i] <- "Australian"
  } else {
    loc[i] <- "Rest of World"
  }
}

loc <- as.data.frame(loc)
colnames(loc) <- c("Location")
rownames(loc) <- sixtyfive_pc$tip.label

# get colour fector for tips
col <- vector()
for (i in 1:length(loc$Location)){
  if(grepl("AUS", loc$Location[i]) || grepl("NZL", loc$Location[i])){
    col[i] <- "#b2182b"
  } else {
    col[i] <- "#4d4d4d"
  }
}

p65 <- ggtree(sixtyfive_pc) +
  geom_cladelabel(node = 58, label = "Rest of \n World", offset = 15, offset.text = 15, fontsize = 5,
                  col = "#b2182b") +
  geom_cladelabel(node = 88, label = "Australian", offset = 15, offset.text = 15, fontsize = 5,
                  col  = "#4d4d4d")

tiff(file = "fig2.tiff", compression = "lzw")
gheatmap(p65, data = loc, width = 0.2, offset = 125,
         colnames = F) +
  scale_fill_manual(values = c("#4d4d4d",  "#b2182b")) +
  theme(legend.text = element_text(size = 12),
  legend.title = element_blank())
dev.off()

## To do:
# consider adding date axis for trees and interpreting time of split