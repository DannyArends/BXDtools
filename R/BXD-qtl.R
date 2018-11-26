#
# BXD-qtl.R
#
# copyright (c) 2017-2020 - Danny Arends, Pjotr Prins, Rob Williams, Gudrun A. Brockmann
# last modified Apr, 2017
# first written Apr, 2017
# 
# Routines to calculate QTL in BXD mice
#

calc.qtl <- function(bxd.genotypes, bxd.phenotypes, pheno.col = 1, n.cores=1, verbose = FALSE) {
  time1 <- proc.time()
  phenotype <- bxd.phenotypes[pheno.col,]
  notNA <- names(which(!is.na(phenotype)))
  phenotype <- phenotype[notNA]
  bxd.subset <- bxd.genotypes[,notNA]
  bxd.numeric <- recode.BXD.genotypes(bxd.subset, rbind(c("D", "1"), c("B", "-1"), c("H", "0"), c("U", "NA")))
  time2 <- proc.time()
  cl <- parallel::makeCluster(n.cores)
  qtls <- parallel::parApply(cl, bxd.numeric, 1, function(genotype) {
    model <- lm(phenotype ~ genotype)
    myAnova <- anova(model)
    return(-log10(myAnova[[5]][1]) * sign(model$coefficients["genotype"]))
  })
  parallel::stopCluster(cl)
  time3 <- proc.time()
  if(verbose){
    cat("Setup:", (time2-time1)[3], "seconds, Mapping:", (time3-time2)[3], "seconds\n")
  }
  invisible(qtls)
}
