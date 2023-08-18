library(ggplot2)
library(reshape2)
library(cowplot)
library(gridGraphics)
library(ape)
library(tidyverse)

## data for decay plot
dat <- read.csv("metadata.csv", header = TRUE)

# get alignments
aln_files <- dir(
  path = "./data/coverage_separated/",
  pattern = ".gatk.fasta$",
  full.names = TRUE
)

aln <- lapply(
  c(aln_files, "./data/COI_653bp_f.fasta"), # NB Inclusion of COI 65% dataset
  function(x) {
    ape::read.dna(x,
      format = "fasta", as.matrix = TRUE, as.character = TRUE
    )
  }
)

names(aln) <- gsub(
  c(aln_files, "COI"),
  pattern = "./data/coverage_separated//|_all.+",
  replacement = ""
)

## Coverage vs time
# cov data
size <- as_tibble(
  read.delim("file_size_filtered.txt", sep = " ", header = FALSE)
  )
colnames(size) <- c("size", "sample")

size <- size %>%
  mutate(megabytes = case_when(
    grepl("G", size) ~ as.numeric(gsub(pattern = "G", replacement = "", size)) * 1000,
    grepl("M", size) ~ as.numeric(gsub(pattern = "M", replacement = "", size)) * 1,
  ))

cov <- left_join(dat, size, by = "sample")


cov_vs_time <- cov %>%
  subset(coverage > 0) %>%
  ggplot() +
  geom_point(
    aes(x = year, y = coverage),
    pch = 21, fill = "dodgerblue", size = 3, alpha = 0.6
  ) +
  ylab("Coverage (%)") +
  xlab("Sample Year") +
  theme_classic(base_size = 13)


cov_vs_size <- cov %>%
  subset(coverage > 0) %>%
  ggplot() +
  geom_point(
    aes(x = megabytes, y = coverage),
    pch = 21, fill = "dodgerblue", size = 3, alpha = 0.6
  ) +
  ylab("") +
  xlab("File Size (Mb)") +
  theme_classic(base_size = 13)

pdf("coverage.pdf", useDingbats = FALSE, width = 8)
cowplot::plot_grid(
  cov_vs_time,
  cov_vs_size,
  labels = "AUTO",
  nrow = 1
)
dev.off()

# genome coverage plots
get_presence_absence <- function(mat) {
  for (i in seq_along(mat[1, ])) {
    for (j in seq_along(mat[, 1])) {
      if (mat[j, i] %in% c("a", "t", "g", "c")) {
        mat[j, i] <- "Base present"
      } else {
        mat[j, i] <- "Base absent"
      }
    }
  }
  print("DOne-1")
  return(melt(mat))
}

pres_abs <- lapply(
  aln,
  function(x) {
    get_presence_absence(x)
  }
)

tile_data <- bind_rows(pres_abs, .id = "id")
colnames(tile_data) <- c("coverage", "sample", "position", "value")

tile_data <- tile_data %>%
  mutate(coverage = gsub(coverage, pattern = "_percent", replacement = ""))

pdf("alignment.pdf", useDingbats = FALSE)
tile_data %>%
  subset(coverage %in% c("5", "25", "50", "65", "COI")) %>%
  ggplot() +
  geom_raster(
    aes(x = position, y = sample, fill = value)
  ) +
  scale_fill_manual(values = c("grey", "dodgerblue"), name = "") +
  scale_y_discrete(labels = NULL, breaks = NULL) +
  facet_wrap(~coverage, scales = "free") +
  theme(
    panel.background = element_rect(fill = "white"),
    legend.position = c(1, 0),
    legend.justification = c(1, 0)
  )
dev.off()
