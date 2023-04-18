# script gets sizes of aus-only clades
library(ape)
library(ggplot2)

source("./analysis_scripts/get_clades.R")

# read in trees
tfiles <- dir(pattern = "./data/relaxed[.]trees") # UPTO: Update tree path
trees <- lapply(tfiles, function(x) read.nexus(x))

# apply burnin
trees <- lapply(trees, function(x) x[-c(1:5000)])

# get aus clade information
# relies on: find.monophyletic(trees[[1]][[1]], tag = "AUS")

# to start, check num aus samples in each
num_aus_total <- c(length(grep(pattern = "AUS", trees[[1]][[1]]$tip.label)),
 length(grep(pattern = "AUS", trees[[2]][[1]]$tip.label)))
num_aus_total
# 221, and 31

max_aus_clades <- list()
for (i in seq_along(trees)) {
max_aus_clade[[i]] <- lapply(trees[[i]],
 function(x) max(lengths(find.monophyletic(x, tag = "AUS"))))
}

max_aus_clade_size <- lapply(max_aus_clade, function(x) unlist(x))
# generating HPD and range
hpd_max_aus_clade_size <- lapply(max_aus_clade_size,
 function(x) quantile(x, probs = c(0.25, 0.975)))
range_max_aus_clade_size <- lapply(max_aus_clade_size, function(x) range(x))

# see node probability before progressing here

# plottingmax aus clade
hist(max_aus_clade / num_aus_total, xlab = "Proportion of Australian Samples",
 main = "Largest Australian-only Clade", col = "dodgerblue")
