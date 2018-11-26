#
# BXD-geneticmap.R
#
# copyright (c) 2017-2020 - Danny Arends, Pjotr Prins, Rob Williams, Gudrun A. Brockmann
# last modified Apr, 2017
# first written Jan, 2017
# 
# Routines to calculate allele frequencies
#

# Calculate the centiMorgan position using Carter-Falconer map locations from the BXD genotype data
calculate.cM.positions <- function(bxd.genotypes, count.Heterozygous = TRUE, start.at.zero = FALSE, verbose = FALSE){
  bxd.map <- attr(bxd.genotypes, "map")
  bxd.map <- cbind(bxd.map, nRecP = 0, nRecS = 0)

  chromosomes <- unique(bxd.map[,"Chr"])
  recombinations <- NULL
  for(chr in chromosomes) {
    onChr <- rownames(bxd.map)[which(bxd.map[,"Chr"] == chr)]
    if(length(onChr) > 1) {                                                               # Make sure there is more then 1 marker on a chromosome
      geno.subset <- bxd.genotypes[onChr, ]                                               # Subset the genotypes on this chromosome
      previousState <- geno.subset[1, ]
      nrecPerChrInd <- rep(0, ncol(bxd.genotypes))
      nrecSinceStart <- 0
      for(x in rownames(geno.subset)[-1]) {                                               # Start at the second marker on the chromosome
        nrecToPrev <- 0
        for(i in 1:ncol(geno.subset)) {
          if(!is.na(geno.subset[x,i]) && !is.na(previousState[i])) {                      # Make sure we don't accidentally compare missing values
            if(previousState[i] != geno.subset[x,i]) {
              if(previousState[i] == "H" || geno.subset[x,i] == "H") {
                if(count.Heterozygous) nrecToPrev <- nrecToPrev + 0.5                     # Count heterozygous genotypes ?
              } else {
                nrecToPrev <- nrecToPrev + 1
                nrecPerChrInd[i] <- nrecPerChrInd[i] + 1
              }
              previousState[i] <- geno.subset[x, i]
            }
          }
          if(is.na(previousState[i]) && !is.na(geno.subset[x, i])) {                      # Is the previous state is an NA, but the current genotype isn't update the previousState
            previousState[i] <- geno.subset[x, i]
          }
        }
        nrecSinceStart = nrecSinceStart + nrecToPrev
        bxd.map[x, "nRecP"] = nrecToPrev
        bxd.map[x, "nRecS"] = nrecSinceStart
        if(verbose) cat(chr, x, nrecToPrev, nrecSinceStart, "\n")
      }
      recombinations <- rbind(recombinations, nrecPerChrInd)
    }
  }
  bxd.map <- cbind(bxd.map, cMraw = NA)                                             # Number of recombinations towards start of chromosome
  bxd.map[,"cMraw"] <- round(100 * (bxd.map[,"nRecS"] / ncol(bxd.map)), digits = 2)        # Remember the RAW cM positions

  R <- (bxd.map[,"nRecP"] / ncol(bxd.genotypes))                                       # Deflate the recombination fractions (KW Broman 11 Aug 2016) 
  r = R/(4 - 6 * R)

  bxd.map <- cbind(bxd.map, imfcf = NA)                                             # Carter-Falconer map distances to previous marker
  bxd.map[,"imfcf"] <- qtl::imf.cf(r)

  bxd.map <- cbind(bxd.map, imfcfsum = NA)                                          # Carter-Falconer map locations to start of chromosome
  for(chr in chromosomes) {
    onChr <- rownames(bxd.map)[which(bxd.map[,"Chr"] == chr)]
    l <- 0
    s <- 0
    if(!start.at.zero) s <- min(as.numeric(bxd.map[onChr, "cM"]))                      # Add the start of the marker on the old genetic map
    l <- l + s
    for(x in 1:length(onChr)) {
      l <- l + bxd.map[onChr[x],"imfcf"]
      bxd.map[onChr[x], "imfcfsum"] <- round(l, digits = 2)
    }
  }
  colnames(recombinations) <- colnames(bxd.genotypes)
  rownames(recombinations) <- chromosomes
  attr(bxd.map, "recombinations") <- recombinations
  return(bxd.map)
}

plot.epochs <- function(bxd.map) {
  recombinations <- attr(bxd.map, "recombinations")
  recombinations <- recombinations[,grep("BXD", colnames(recombinations))]
  chromosomes <- rownames(recombinations)
  op <- par(mar = c(5, 4, 1, 1))
  cols <- colorRampPalette(RColorBrewer::brewer.pal(6, "Dark2"))(nrow(recombinations))

  cohorts <- c(which(colnames(recombinations) == "BXD42")+0.5,
               which(colnames(recombinations) == "BXD102")+0.5,
               which(colnames(recombinations) == "BXD157")+0.5,
               which(colnames(recombinations) == "BXD186")+0.5, ncol(recombinations)+0.5)
  pY <- apply(recombinations, 2, sum)

  plot(x = c(0.5, 0.5 + ncol(recombinations)), y = c(0, max(pY, na.rm = TRUE) * 1.2), t = 'n', xlab = "BXD Epoch", ylab = "# Recombinations", 
       xaxt='n', xaxs="i", yaxs="i", las=2, cex.lab=1.5,cex.axis=1.5)

  for(x in nrow(recombinations):1){
    rect(rep(1:ncol(recombinations)) - 0.5, pY, rep(1:ncol(recombinations)) + 0.5, pY - recombinations[x,], col=cols[x], border = "white", lwd =0.1)
    pY <- pY - recombinations[x,]
  }
  abline(v = cohorts)
  legend_order <- matrix(1:20,ncol=2, byrow = FALSE)
  axis(1, at = midpoints(cohorts), c("Taylor '70 - '90", "Williams '90", "Williams '08", "Williams & Palmer '10", "Williams '14"),cex.axis=1.5)
  box()
  legend("topleft", paste0("Chr ", chromosomes)[legend_order], fill=cols[legend_order], cex=1, ncol=2, border = "white")
}

plot.map <- function(bxd.genotypes, add.markers = TRUE, highlight.markers = c(), main = "") {
  geneticMap <- attr(bxd.genotypes, "map")
  chromosomes <- unique(geneticMap[,"Chr"])
  maxY <- 1.15 * max(as.numeric(unique(geneticMap[,"Mb"])))
  plot(x = c(1, length(chromosomes)), y = c(0, maxY), t = 'n', xaxt='n', yaxt = 'n', xlab="Chromosome", ylab="Position (Mb)", main = main)
  axis(1, at = 1:length(chromosomes), chromosomes)
  axis(2, at = seq(0, maxY, 10), seq(0, maxY, 10), las = 2)
  
  for(chr.n in 1:length(chromosomes)){
    chr.map <- geneticMap[geneticMap[,"Chr"] == chromosomes[chr.n], ]
    chr.length <- max(as.numeric(chr.map[, "Mb"]))
    segments(chr.n, 0, chr.n, chr.length)
    if(add.markers) {
      colz <- rep("black", nrow(chr.map))
      colz[which(rownames(chr.map) %in% highlight.markers)] <- "green"
      points(x = rep(chr.n, nrow(chr.map)), y = as.numeric(chr.map[,"Mb"]), pch = 95, col=colz)
    }
  }
}

plot.qtl <- function(bxd.genotypes, lodscores, chr = NULL, highlight.markers = c(), gap = 25, main = "", type = 'l', ...) {
  geneticMap <- attr(bxd.genotypes, "map")
  chromosomes <- unique(geneticMap[,"Chr"])
  if(!is.null(chr)){
    chromosomes <- chromosomes[which(chromosomes %in% chr)]
  }
  total.chr.length <- -gap
  chr.start.pos <- c()
  chr.lengths <- c()
  
  # Fix when lodscores do not have names()  
  if(is.null(names(lodscores)) && length(lodscores) == nrow(geneticMap)) {
    names(lodscores) <- rownames(geneticMap)
  }
  
  for(chr.n in 1:length(chromosomes)){
    total.chr.length <- total.chr.length + gap
    chr.start.pos <- c(chr.start.pos, total.chr.length)
    chr.map <- geneticMap[geneticMap[,"Chr"] == chromosomes[chr.n], ]
    chr.length <- max(as.numeric(chr.map[, "Mb"]))
    chr.lengths <- c(chr.lengths, chr.length)
    total.chr.length <- total.chr.length + chr.length
  }
  plot(x=c(0, total.chr.length), y = c(min(lodscores,na.rm = TRUE), max(lodscores,na.rm = TRUE)), t = 'n', xlab="Position", ylab="-log10(pvalue)", main = main, xaxt='n')
  for(chr.n in 1:length(chromosomes)){
    chr.map <- geneticMap[geneticMap[,"Chr"] == chromosomes[chr.n], ]
    points(x=as.numeric(chr.map[, "Mb"]) + chr.start.pos[chr.n], y = rep(0, nrow(chr.map)), pch = "|", cex = 0.3)
    colz <- rep("black", nrow(chr.map))
    colz[which(rownames(chr.map) %in% highlight.markers)] <- "green"
    points(x = as.numeric(chr.map[, "Mb"]) + chr.start.pos[chr.n], y = lodscores[rownames(chr.map)], type = type, ...)
  }
  axis(1, at = (chr.lengths / 2.0) + chr.start.pos, chromosomes)
  abline(h = -log10(0.05 / length(lodscores)), lty=2, col="orange")
  abline(h = -log10(0.01 / length(lodscores)), lty=2, col="green")
}
