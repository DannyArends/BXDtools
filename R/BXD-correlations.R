#
# BXD-correlations.R
#
# copyright (c) 2017-2020 - Danny Arends, Pjotr Prins, Rob Williams, Gudrun A. Brockmann
# last modified Apr, 2017
# first written Apr, 2017
# 
# Routines to calculate correlations between phenotypes in BXD mice
#

calc.correlations <- function(phenotype, bxd.phenotypes, top = 100, matched = 0.9, method = c("pearson", "kendall", "spearman")) {
  bxd.phenotypes <- phenotype.add.class(bxd.phenotypes)
  phenotype.descriptions <- attr(bxd.phenotypes, "annotation")
  bxd.subset <- bxd.phenotypes[,names(phenotype)]
  perc.missing <- apply(bxd.subset,1,function(x){
    return(sum(is.na(x))/length(x))
  })
  bxd.subset <- bxd.subset[-which(perc.missing >= matched),]
  correlations <- cor(phenotype, t(bxd.subset), use = "pairwise.complete.obs", method = method[1])
  names(correlations) <- rownames(bxd.subset)
  sorted <- sort(abs(correlations), decreasing = TRUE)[1:top]
  return(data.frame(cbind( phenotype.descriptions[names(sorted),], "R" = round(correlations[names(sorted)], 3))))
}

test.correlation <- function(){
  require("BXDtools")

  data('bxd.phenotypes', package='BXDtools', envir=environment())
  bxd.phenotypes <- get("bxd.phenotypes", envir=environment()) 

  bxd.phenotypes <- phenotype.add.class(bxd.phenotypes)

  perc.missing <- apply(bxd.phenotypes, 2, function(x){
    return(sum(is.na(x))/length(x))
  })

  N <- 30
  availableBXD <- names(which(perc.missing < 0.9))
  phenotype <- rnorm(N)
  names(phenotype) <- sample(availableBXD, N)

  calc.correlations(phenotype, bxd.phenotypes, 5)
  
  table(attr(bxd.phenotypes, "annotation")[,"class"])
}
