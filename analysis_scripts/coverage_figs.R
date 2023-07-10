library(ggplot2)
library(reshape2)
library(cowplot)
library(gridGraphics)
library(ape)
library(tidyverse)

## data for decay plot
dat <- read.csv("./data/tableS1.csv", header = TRUE)

## Doing Coverage plots
# COI coverage
mat <- read.dna(
  file = "./data/COI_653bp_f.fasta",
  format = "fasta", as.matrix = TRUE, as.character = TRUE
)

for (i in seq_along(mat[1, ])) {
  for (j in seq_along(mat[, 1])) {
    if (mat[j, i] %in% c("a", "t", "g", "c")) {
      mat[j, i] <- "Base present"
    } else {
      mat[j, i] <- "Base absent"
    }
  }
}
cov_coi <- melt(mat)

# 65pc coverage
mat <- read.dna(
  file = "./data/65percent_mitogenomes.fasta",
  format = "fasta",
  as.matrix = TRUE,
  as.character = TRUE)

for (i in seq_along(mat[1, ])) {
  for (j in seq_along(mat[, 1])) {
    if (mat[j, i] %in% c("a", "t", "g", "c")) {
      mat[j, i] <- "Base present"
    } else {
      mat[j, i] <- "Base absent"
    }
  }
}
cov_65 <- melt(mat)

# 5pc coverage
mat <- read.dna(
  file = "./data/5percent_mitogenomes.fasta",
  format = "fasta",
  as.matrix = TRUE,
  as.character = TRUE
)

for (i in seq_along(mat[1, ])) {
  for (j in seq_along(mat[, 1])) {
    if (mat[j, i] %in% c("a", "t", "g", "c")) {
      mat[j, i] <- "Base present"
    } else {
      mat[j, i] <- "Base absent"
    }
  }
}
cov_5 <- melt(mat)

# Panel figure with all of these

decay <- ggplot(dat) +
  geom_point(aes(x = Year.of.collection, y = Percent.Mitogenome.Coverage),
    pch = 21, bg = alpha("steelblue", 0.5), size = 3
  ) +
  xlab("Year of collection") +
  ylab("Mitogenome coverage (%)") +
  ylim(0, 100) +
  theme_classic(base_size = 13)

pc5 <- ggplot(cov_5, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  xlab("Mitogenome Postion") +
  ylab("Individual") +
  scale_fill_manual(values = c("#f0f0f0", "#636363"), name = "") +
  scale_y_discrete(labels = NULL, breaks = NULL) +
  theme(
    panel.background = element_rect(fill = "white"),
    legend.position = "none"
  )

pc65 <- ggplot(cov_65, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  xlab("Mitogenome Postion") +
  ylab("Individual") +
  scale_fill_manual(values = c("#f0f0f0", "#636363"), name = "") +
  scale_y_discrete(labels = NULL, breaks = NULL) +
  theme(
    panel.background = element_rect(fill = "white"),
    legend.position = "none"
  )

legend <- get_legend(
  ggplot(cov_coi, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  xlab("COI Postion") +
  ylab("") +
  scale_fill_manual(values = c("#f0f0f0", "#636363"), name = "") +
  scale_y_discrete(labels = NULL, breaks = NULL) +
  theme(
    panel.background = element_rect(fill = "white"),
    legend.position = "none"
  )
)
coi <- ggplot(cov_coi, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill = value)) +
  xlab("COI Postion") +
  ylab("") +
  scale_fill_manual(values = c("#f0f0f0", "#636363"), name = "") +
  scale_y_discrete(labels = NULL, breaks = NULL) +
  theme(
    panel.background = element_rect(fill = "white"),
    legend.position = "none"
  )


panel <- plot_grid(
  decay, 
  pc5, 
  pc65, 
  coi, 
  labels = c("a", "b", "c", "d"), 
  nrow = 2, 
  ncol = 2)

plot_grid(panel, legend, rel_widths = c(2, .4))


## file size  and coverge plot
size <- as_tibble(read.delim("file_size_filtered.txt", sep = " ", header = FALSE))
colnames(size) <- c("size", "sample")

cov <- as_tibble(read.delim("coverage.txt", sep = " ", header = FALSE))
colnames(cov) <- c("sample", "coverage")

# some mismatch between samples. right_join ensures all file sizes defined
data <- left_join(cov, size, by = "sample") %>%
  mutate(megabytes = case_when(
    grepl("G", size) ~ as.numeric(gsub(pattern = "G", replacement = "", size)) * 1000,
    grepl("M", size) ~ as.numeric(gsub(pattern = "M", replacement = "", size)) * 1,
  ))

# cov vs. size
svg("coverage_vs_filesize.svg")
  ggplot(data) +
    geom_point(
      aes(x = megabytes, y = 100 * coverage),
      pch = 21,
      fill = "steelblue",
      size = 3,
      alpha = 0.6
    ) +
    ylim(0, 100) +
    ylab("% Mitogenome Coverage") +
    theme_minimal()
dev.off()

# coverage threshold plot

totalSamples <- length(data$sample)
all <- data %>%
  select(coverage) %>%
  group_by(coverage) %>%
  summarise(n = n()) %>%
  mutate(n = totalSamples - cumsum(n))

ausSamples <- length(which(grepl("AM|Au_", data$sample)))
aus <- data %>%
  filter(grepl("AM|Au_", sample)) %>%
  select(coverage) %>%
  group_by(coverage) %>%
  summarise(n = n()) %>%
  mutate(n = ausSamples - cumsum(n))



svg("number_above_coverage.svg")
bind_rows("Aus Samples" = aus, "All Samples" = all, .id = "group") %>%
  ggplot() +
  geom_line(
    aes(x = 100 * coverage, y = n, col = group),
  ) +
  scale_color_manual(values = c("black", "darkred")) +
  xlab("Mitogenome Coverage Threshold (%)") +
  ylab("Number of Samples above Threshold") +
  annotate(
    "segment",
    x = c(5, 50, 95),
    xend = c(5, 50, 95),
    y = c(0, 0, 0),
    yend = c(Inf, Inf, Inf),
    colour = "blue"
  ) +
  annotate(
    "label",
    x = c(5, 50, 95),
    y = c(300, 300, 300),
    label = c("5%", "50%", "95%")
  ) +
  theme_minimal()
dev.off()



