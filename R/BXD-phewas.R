#
# BXD-phewas.R
#
# copyright (c) 2017-2020 - Danny Arends, Pjotr Prins, Rob Williams, Gudrun A. Brockmann
# last modified Apr, 2017
# first written Jan, 2017
# 
# Routines to perform pheWAS on the BXD dataset
#

do.BXD.phewas <- function(bxd.genotypes, bxd.phenosomes, marker = "rs31443144", LRS = FALSE) {
  scores <- apply(bxd.phenosomes, 1, function(pheno, geno) {
    result <- tryCatch({
      myanova <- anova(lm(pheno ~ geno))
      if(!LRS) return(-log10(myanova[[5]][1]))
      return(sum(!is.na(pheno)) * log(sum(myanova$"Sum Sq") / myanova$"Sum Sq"[2]))
    }, error = function(e) {
      cat("Error caught at marker: '", marker, "'\n")
      return(NA)
    })
    return(result)
  }, bxd.genotypes[marker,])
  attr(scores, "marker") <- marker
  attr(scores, "LRS") <- LRS
  return(scores)
}

do.BXD.phewas.all <- function(bxd.genotypes, bxd.phenosomes) {
  res <- lapply(rownames(bxd.genotypes), function(marker.name) {
    cat("marker: ", marker.name, "\n")
    do.BXD.phewas(bxd.genotypes, bxd.phenosomes, marker.name)
  })
  return(res)
}

plot.phewas <- function(pvalues, bxd.phenosomes, do.sort = FALSE, decreasing = FALSE,
                        main = "BxD PheWAS results", pch = 19, cex = 0.6, significance = 19, type = "h",
                        colorSeed = 1, colorRange = c("darkslateblue", "hotpink1", "forestgreen", "orange", "black", "firebrick1"), ...) {
  classes <- attr(bxd.phenosomes, "annotation")[,"class"]
  marker <- attr(scores, "marker")
  ylab <- "-log10(P)"
  if(attr(scores, "LRS")) ylab = "LRS"
  if(do.sort) {
    ordering <- c()
    for(pheclass in unique(classes)){
      indexes <- which(classes == pheclass)
      # Use Radix sorting and put the NA scores last
      sorting <- sort(scores[indexes], decreasing = decreasing, index.return = TRUE, method = "quick", na.last=NA)
      ordering <- c(ordering, names(sorting$x))
    }
    scores <- scores[ordering]
    classes.inOrder <- classes[ordering]
  }else{
    classes.inOrder <- classes
  }
  # Set the colorSeed to NA, to cycle colors to groups each plot
  if(!is.na(colorSeed)) set.seed(colorSeed)
  # Generate colors from colorRange and mix them up for beter visability
  mcolors <- sample(colorRampPalette(colorRange)(length(unique(classes))))
  names(mcolors) <- unique(classes)

  plot(scores, col = mcolors[classes.inOrder], pch = pch, 
                        cex = cex, type = type, xaxt='n', xlab = "Phenosome", ylab = ylab, main=main, ...)
  legend("topright", paste0(names(mcolors), " (N = ", table(classes)[names(mcolors)], ")"), 
                     col = mcolors, pch = pch, cex = cex, ncol = 2)
  abline(h = significance, lty=2)
  results <- cbind(classes.inOrder, scores)
  significant <- results[which(as.numeric(results[,2]) > significance),]
  legend("topleft", paste0(rownames(significant), " ", significant[,1], " ", ylab, "=", round(as.numeric(significant[,2]),1)), cex = cex, fill = mcolors[significant[,1]])

  invisible(results)
}
