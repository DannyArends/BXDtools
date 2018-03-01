#
# BXD-phewas.R
#
# copyright (c) 2017-2020 - Danny Arends, Pjotr Prins, Rob Williams, Gudrun A. Brockmann
# last modified Apr, 2017
# first written Jan, 2017
# 
# Routines to perform pheWAS on the BXD dataset
#

do.BXD.phewas <- function(bxd.genotypes, bxd.phenosomes, marker = "rs31443144") {
  pvalues <- apply(bxd.phenosomes, 1, function(pheno, geno) {
    result <- tryCatch({
      anova(lm(pheno ~ geno))[[5]][1]
    }, error = function(e) {
      cat("Error caught at marker: '", marker, "'\n")
      return(NA)
    })
    return(result)
  }, bxd.genotypes[marker,])
  attr(pvalues, "marker") <- marker
  return(pvalues)
}

do.BXD.phewas.all <- function(bxd.genotypes, bxd.phenosomes) {
  res <- lapply(rownames(bxd.genotypes), function(marker.name) {
    cat("marker: ", marker.name, "\n")
    do.BXD.phewas(bxd.genotypes, bxd.phenosomes, marker.name)
  })
  return(res)
}

plot.phewas <- function(pvalues, bxd.phenosomes, do.sort = FALSE, decreasing = FALSE,
                        main = "BxD PheWAS results", pch = 19, cex = 0.6, type = "h",
                        colorSeed = 1, colorRange = c("darkslateblue", "hotpink1", "forestgreen", 
                                                      "orange", "black", "firebrick1")) {
  classes <- attr(bxd.phenosomes, "annotation")[,"class"]
  marker <- attr(pvalues, "marker")
  if(do.sort) {
    ordering <- c()
    for(pheclass in unique(classes)){
      indexes <- which(classes == pheclass)
      # Use Radix sorting and put the NA pvalues last
      sorting <- sort(pvalues[indexes], decreasing = decreasing, index.return = TRUE, method = "quick", na.last=NA)
      ordering <- c(ordering, indexes[sorting$ix])
    }
    pvalues <- pvalues[ordering]
    classes.inOrder <- classes[ordering]
  }else{
    classes.inOrder <- classes
  }
  # Set the colorSeed to NA, to cycle colors to groups each plot
  if(!is.na(colorSeed)) set.seed(colorSeed)
  # Generate colors from colorRange and mix them up for beter visability
  mcolors <- sample(colorRampPalette(colorRange)(length(unique(classes))))
  names(mcolors) <- unique(classes)

  main <- paste0("Phewas on ", length(classes), " phenotypes, at ", marker)
  plot(-log10(pvalues), col = mcolors[classes.inOrder], pch = pch, 
                        cex = cex, type = type, xaxt='n', xlab = "Phenosome", ylab="-log10(P)", main=main)
  legend("topright", paste0(names(mcolors), " (N = ", table(classes)[names(mcolors)], ")"), 
                     col = mcolors, pch = pch, cex = cex)
  invisible(cbind(classes.inOrder, pvalues))
}

