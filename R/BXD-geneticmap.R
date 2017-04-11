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
  geneticMap <- attr(bxd.genotypes, "map")
  geneticMap <- cbind(geneticMap, nRecP = 0, nRecS = 0)

  chromosomes <- unique(geneticMap[,"Chr"])
  for(chr in chromosomes) {
    onChr <- rownames(geneticMap)[which(geneticMap[,"Chr"] == chr)]
    if(length(onChr) > 1) {                                                               # Make sure there is more then 1 marker on a chromosome
      geno.subset <- bxd.genotypes[onChr, ]                                               # Subset the genotypes on this chromosome
      previousState <- geno.subset[1, ]
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
              }
              previousState[i] <- geno.subset[x, i]
            }
          }
          if(is.na(previousState[i]) && !is.na(geno.subset[x, i])) {                      # Is the previous state is an NA, but the current genotype isn't update the previousState
            previousState[i] <- geno.subset[x, i]
          }
        }
        nrecSinceStart = nrecSinceStart + nrecToPrev
        geneticMap[x, "nRecP"] = nrecToPrev
        geneticMap[x, "nRecS"] = nrecSinceStart
        if(verbose) cat(chr, x, nrecToPrev, nrecSinceStart, "\n")
      }
    }
  }
  geneticMap <- cbind(geneticMap, cMraw = NA)                                             # Number of recombinations towards start of chromosome
  geneticMap[,"cMraw"] <- round(100 * (geneticMap[,"nRecS"] / ncol(geneticMap)),2)        # Remember the RAW cM positions

  R <- (geneticMap[,"nRecP"] / ncol(bxd.genotypes))                                       # Deflate the recombination fractions (KW Broman 11 Aug 2016) 
  r = R/(4 - 6 * R)

  geneticMap <- cbind(geneticMap, imfcf = NA)                                             # Carter-Falconer map distances to previous marker
  geneticMap[,"imfcf"] <- qtl::imf.cf(r)

  geneticMap <- cbind(geneticMap, imfcfsum = NA)                                          # Carter-Falconer map locations to start of chromosome
  for(chr in chromosomes) {
    onChr <- rownames(geneticMap)[which(geneticMap[,"Chr"] == chr)]
    l <- 0
    s <- 0
    if(!start.at.zero) s <- min(as.numeric(geneticMap[onChr, "cM"]))                      # Add the start of the marker on the old genetic map
    l <- l + s
    for(x in 1:length(onChr)) {
      l <- l + geneticMap[onChr[x],"imfcf"]
      geneticMap[onChr[x], "imfcfsum"] <- round(l, d = 2)
    }
  }
  return(geneticMap)
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
      idx2 <- which(colz == "green")
      points(x = rep(chr.n, length(idx2)), y = as.numeric(chr.map[idx2,"Mb"]), pch = 95, col="green")
    }
  }
}

plot.qtl <- function(bxd.genotypes, lodscores, highlight.markers = c(), gap = 25, main = "") {
  geneticMap <- attr(bxd.genotypes, "map")
  chromosomes <- unique(geneticMap[,"Chr"])
  total.chr.length <- -gap
  chr.start.pos <- c()
  
  # Fix when lodscores do not have names()  
  if(is.null(names(lodscores)) && length(lodscores) == nrow(geneticMap)) {
    names(lodscores) <- rownames(geneticMap)
  }
  
  for(chr.n in 1:length(chromosomes)){
    total.chr.length <- total.chr.length + gap
    chr.start.pos <- c(chr.start.pos, total.chr.length)
    chr.map <- geneticMap[geneticMap[,"Chr"] == chromosomes[chr.n], ]
    chr.length <- max(as.numeric(chr.map[, "Mb"]))
    total.chr.length <- total.chr.length + chr.length
  }
  plot(x=c(0, total.chr.length), y = c(min(lodscores,na.rm = TRUE), max(lodscores,na.rm = TRUE)), t = 'n', xlab="Cumulative position (Mb)", ylab="-log10(pvalue)", main = main)
  for(chr.n in 1:length(chromosomes)){
    chr.map <- geneticMap[geneticMap[,"Chr"] == chromosomes[chr.n], ]
    points(x=as.numeric(chr.map[, "Mb"]) + chr.start.pos[chr.n], y = rep(0, nrow(chr.map)), pch = "|", cex = 0.3)
    colz <- rep("black", nrow(chr.map))
    colz[which(rownames(chr.map) %in% highlight.markers)] <- "green"
    points(x=as.numeric(chr.map[, "Mb"]) + chr.start.pos[chr.n], y = lodscores[rownames(chr.map)], t = 'l')
    points(x=as.numeric(chr.map[, "Mb"]) + chr.start.pos[chr.n], y = lodscores[rownames(chr.map)], col=colz, t = 'p',pch = 18)
  }
  abline(h = -log10(0.05 / length(lodscores)), lty=2, col="orange")
  abline(h = -log10(0.01 / length(lodscores)), lty=2, col="green")
}
