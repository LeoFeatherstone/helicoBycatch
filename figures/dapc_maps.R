## Mapping DAPC output
library(oz)
library(ggplot2)

setwd(paste0(getwd(), "/figures"))

# 5pc map
five_mito_problocations <- read.csv("5pc_mito_locations_probs.csv")

post <- five_mito_problocations[,2:3]
post <- post[order(post[,1]),]

col <- vector()
for (i in 1:length(five_mito_problocations$Subspecies)){
  if (five_mito_problocations$Subspecies[i] == "H. a. armigera"){
    col[i] <- "#b2182b"
  } else {
    col[i] <- "#4d4d4d"
  }
}

# Making barplot 
tiff(file = "figS1A.tiff", compression = "lzw")
  layout(mat = matrix(c(1, 3, 2, 3), 
                      nrow = 2, 
                      ncol = c(2,1)), 
        heights = c(1, 2),        # Heights of the two rows
        widths = c(1.5, 1.5))     # Widths of the two columns

  plot.new()
  legend('center',
        legend = c("Australia ",
                    "Rest of World "),
        col = alpha(col, 0.55),
        pch = 16,
        cex = 1.5,
        pt.cex = 2.5,
        border = F)

  post <- five_mito_problocations[1:221,2:3]
  post <- post[,c(2,1)]
  bar <- barplot(t(post[order(post[,2]),]), 
                col = alpha(col, 0.55),
                border= NA,
                space = 0,
                xlab="Individuals", 
                ylab="Posterior Probability",
                names.arg = rep(" ", length(post[,1])),
                font.lab=2,
                las = 1)

  # Making Map
  oz <- oz(states = T)
  points(jitter(five_mito_problocations$Lon_final, amount = 0.6), 
        jitter(five_mito_problocations$Lat_final, 
        amount = 0.6),
        cex = 2.5, pch = 16, col = alpha(col, 0.55))
dev.off()

# Mapping COI Data
locsnew2 <- read.csv("COI653_locs_TEMP.csv")
post <- locsnew2[,3:4]
post <- post[order(post[,1]),]

col <- vector()
for (i in 1:length(locsnew2$Subspecies)){
  if (locsnew2$Subspecies[i] == "H. a. armigera"){
    col[i] <- "#b2182b"
  } else {
    col[i] <- "#4d4d4d"
  }
}

# Making barplot 
tiff(file = "figS1B.tiff", compression = "lzw")
  layout(mat = matrix(c(1, 3, 2, 3), 
                      nrow = 2, 
                      ncol = c(2,1)),
        heights = c(1, 2),    # Heights of the two rows
        widths = c(1.5, 1.5))     # Widths of the two columns

  plot.new()
  legend('center',
                legend = c("Australia ",
                          "Rest of World "),
                col = alpha(unique(col), 0.55),
                pch = 16,
                cex = 1.5,
                pt.cex = 2.5,
                border = F)


  bar <- barplot(t(post), 
          col = alpha(unique(col), 0.55),
          border= NA,
          space = 0,
          xlab="Individuals", 
          ylab="Posterior Probability",
          names.arg = rep(" ", length(post[,1])),
          font.lab=2,
          las = 1)

  # Making Map
  oz <- oz(states = T)
  oz <- points(jitter(locsnew2$Longitude, amount = 0.6), jitter(locsnew2$Latitude, 
                                                          amount = 0.6),
        cex = 2.5, pch = 16, col = alpha(col, 0.55))
dev.off()


# Mapping 65% Data for preprint
dapc_65_locs <- read.csv("65pc_dapc_locations.csv")

post <- dapc_65_locs[,c(4,5)]
post <- post[order(post[,1]),]

  col <- vector()
  for (i in 1:length(dapc_65_locs$Subspecies)){
    if (dapc_65_locs$Subspecies[i] == "H. a. armigera"){
      col[i] <- "#b2182b"
    } else {
      col[i] <- "#4d4d4d"
    }
  }

  # Making barplot 
tiff(file = "fig2.tiff", compression = "lzw")
  layout(mat = matrix(c(1, 3, 2, 3), 
                      nrow = 2, 
                      ncol = c(2,1)), 
        heights = c(1, 2),        # Heights of the two rows
        widths = c(1.5, 1.5))     # Widths of the two columns

  plot.new()
  legend('center',
        legend = c("Australia", "Rest of World"),
        col = alpha(c("#4d4d4d", "#b2182b"), 0.55),
        pch = 16,
        cex = 1.5,
        pt.cex = 2.5,
        border = F)

  bar <- barplot(t(post), 
                col = alpha(col, 0.55),
                border= NA,
                space = 0,
                xlab="Individuals", 
                ylab="Posterior Probability",
                names.arg = rep(" ", length(post[,1])),
                font.lab=2,
                las = 1)

  # Making Map
  oz <- oz(states = T)
    points(jitter(dapc_65_locs$Longitude, amount = 0.6), 
        jitter(dapc_65_locs$Latitude, 
                amount = 0.6),
        cex = 2.5, pch = 16, col = alpha(col, 0.55))
dev.off()

















