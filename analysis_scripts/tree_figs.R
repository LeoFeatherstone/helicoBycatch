## Plotting for figures 2 and sups
library(ggtree)
library(ggplot2)
library(treeio)

# 5pc
five_pc <- treeio::read.beast("./data/5percent_mitogenomes_relaxed_MCC.trees")

# get locations for tips
loc <- vector()
for (i in seq_along(five_pc@phylo$tip.label)) {
  if (
    grepl("AUS_", five_pc@phylo$tip.label[i])
    ||
    grepl("NZL_", five_pc@phylo$tip.label[i])
  ) {
    loc[i] <- "Conferta"
  } else {
    loc[i] <- "Armigera"
  }
}

loc <- as.data.frame(loc)
colnames(loc) <- c("Location")
rownames(loc) <- five_pc@phylo$tip.label

# get colour vector for tips
col <- vector()
for (i in seq_along(loc$Location)) {
  if (grepl("AUS", loc$Location[i]) || grepl("NZL", loc$Location[i])) {
    col[i] <- "#b2182b"
  } else {
    col[i] <- "#4d4d4d"
  }
}

# setting up basic tree plot
p5 <- ggtree(five_pc, aes(color = posterior)) +
  scale_color_continuous(name = "Posterior\nProbability") +
  geom_treescale(x = 0, y = 100, label = "yrs") +
  theme(legend.position = "left")

gheatmap(p5,
  data = loc, width = 0.2, offset = 650,
  colnames = F
) +
  scale_fill_manual(values = c("#b2182b", "#4d4d4d"), name = "") +
  theme(legend.text = element_text(size = 12)) +
  scale_x_ggtree(breaks = c(0, 1000, 2000))

# 65 pc  # have to clear history in all windows before this will work.
sixtyfive_pc <- treeio::read.beast(
  "./data/65percent_mitogenomes_relaxed_MCC.trees"
)

# get locations for tips
loc <- vector()
for (i in seq_along(sixtyfive_pc@phylo$tip.label)) {
  if (
    grepl("AUS_", sixtyfive_pc@phylo$tip.label[i])
    ||
    grepl("NZL_", sixtyfive_pc@phylo$tip.label[i])
  ) {
    loc[i] <- "Conferta"
  } else {
    loc[i] <- "Armigera"
  }
}

loc <- as.data.frame(loc)
colnames(loc) <- c("Location")
rownames(loc) <- sixtyfive_pc@phylo$tip.label

# get colour vector for tips
col <- vector()
for (i in seq_along(loc$Location)) {
  if (grepl("AUS", loc$Location[i]) || grepl("NZL", loc$Location[i])) {
    col[i] <- "#b2182b"
  } else {
    col[i] <- "#4d4d4d"
  }
}

# setting up basic tree plot
p65 <- ggtree(sixtyfive_pc, aes(color = posterior)) +
  scale_color_continuous(name = "Posterior\nProbability") +
  geom_treescale(x = 0, y = 50, label = "yrs") +
  theme(legend.position = "left")

gheatmap(p65,
  data = loc, width = 0.2, offset = 50,
  colnames = FALSE
) +
  scale_fill_manual(values = c("#b2182b", "#4d4d4d"), name = "") +
  theme(legend.text = element_text(size = 12)) +
  scale_x_ggtree(breaks = seq(from = 1600, to = 2000, by = 100))

############## COI Tree plot below for fig 4 ###################

coi_tree <- read.newick("./data/COI_653bp_f.treefile")

loc <- vector()
for (i in seq_along(coi_tree$tip.label)){
  if (
    grepl("AUS_", coi_tree$tip.label[i])
    ||
    grepl("NZL_", coi_tree$tip.label[i])
  ){
    loc[i] <- "Conferta"
  } else {
    loc[i] <- "Armigera"
  }
}

# get locations for tips

loc <- as.data.frame(loc)
colnames(loc) <- c("Location")
rownames(loc) <- coi_tree$tip.label

# get colour vector for tips
col <- vector()
for (i in seq_along(loc$Location)){
  if(grepl("Conferta", loc$Location[i])) {
    col[i] <- "#b2182b"
  } else {
    col[i] <- "#4d4d4d"
  }
}

g <- factor(loc$Location)

group_info <- split(coi_tree$tip.label, g)
test <- groupOTU(coi_tree, group_info)
ggtree(test, aes(color = group), branch.length = "none", layout = "circular") +
  scale_color_manual(values  = c("#b2182b", "#4d4d4d"), name = "") +
  theme(legend.text = element_text(size  = 12)) +
  geom_tiplab(size = 1)
