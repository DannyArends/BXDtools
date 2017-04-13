#
# BXD-correlations.R
#
# copyright (c) 2017-2020 - Danny Arends, Pjotr Prins, Rob Williams, Gudrun A. Brockmann
# last modified Apr, 2017
# first written Apr, 2017
# 
# Routines to calculate correlations between phenotypes in BXD mice
#

calc.correlations <- function(phenotype, bxd.phenotypes, top = 100, method = c("pearson", "kendall", "spearman")) {
  bxd.subset <- bxd.phenotypes[,names(phenotype)]
  correlations <- cor(phenotype, t(bxd.subset), use = "pairwise.complete.obs", method = method[1])
  names(correlations) <- rownames(bxd.subset)
  sorted <- sort(abs(correlations), decreasing = TRUE)[1:top]
  cbind(names(sorted), correlations[names(sorted)])
}

test.correlation <- function(){
  availableBXD <- colnames(bxd.phenotypes)[3:100]
  phenotype <- rnorm(30)
  names(phenotype) <- sample(availableBXD, 30)

}

