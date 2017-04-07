

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
                        colorSeed = 1, colorRange = c("darkslateblue", "hotpink1", "forestgreen", "orange", "firebrick1")) {
  classes <- attr(bxd.phenosomes, "annotation")[,"class"]
  if(do.sort) {
    ordering <- c()
    for(pheclass in unique(classes)){
      indexes <- which(classes == pheclass)
      # Use Radix sorting and put the NA pvalues last
      sorting <- sort(pvalues[indexes], decreasing = decreasing, index.return = TRUE, method = "radix", na.last=TRUE)
      ordering <- c(ordering, indexes[sorting$ix])
    }
    pvalues <- pvalues[ordering]
    classes.inOrder <- classes[ordering]
  }else{
    classes.inOrder <- classes
  }
  if(!is.na(colorSeed)) set.seed(colorSeed)                                   # Set the colorSeed to NA, to cycle colors to groups each plot
  mcolors <- sample(colorRampPalette(colorRange)(length(unique(classes))))    # Generate colors from colorRange and mix them up for beter visability
  names(mcolors) <- unique(classes)

  plot(-log10(pvalues), col = mcolors[classes.inOrder], pch = pch, cex = cex, type = type, xaxt='n', xlab = "Phenosome", ylab="-log10(P)")
  legend("topright", names(mcolors), col = mcolors, pch = pch, cex = cex)
  return(cbind(classes.inOrder, pvalues))
}