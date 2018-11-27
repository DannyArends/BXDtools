#
# BXD-utils.R
#
# copyright (c) 2017-2020 - Danny Arends, Pjotr Prins, Rob Williams, Gudrun A. Brockmann
# last modified Apr, 2017
# first written Jan, 2017
# 
# Additional functions for the BXDtools package
#

midpoints <- function(x, s = 1){
  p <- NULL
  for(i in 1:length(x)){
    p <- c(p, (s + x[i] / 2))
    s <- (p[i] - s)
  }
  return(p)
}

test.BXD <- function(x) {
  require("BXDtools")
  #bxd.geno <- download.BXD.genotypes()
  #bxd.genotypes <- only.BXD.genotypes(bxd.geno)
  #bxd.genotypes <- add.PandF1(bxd.genotypes)

  data('bxd.genotypes', package='BXDtools', envir=environment())
  bxd.genotypes <- get("bxd.genotypes", envir=environment()) 
  bxd.map <- calculate.cM.positions(bxd.genotypes)

  plot.epochs(bxd.map)

  calc.alleleFrequency(bxd.genotypes[1,])

  pM <- NULL
  for(x in rownames(attr(bxd.genotypes,"map"))[1:10]){
    pvals <- c()
    for(y in rownames(attr(bxd.genotypes,"map"))[1:1000]){
      pvals <- c(pvals, calc.LD(bxd.genotypes, mname1 = x, mname2 = y)[["P-value"]])
    }
    pM <- rbind(pM, pvals)
  }
  
  bxd.maf.genotypes <- recode.BXD.to.MAF(bxd.genotypes)
  bxd.numeric <- recode.BXD.genotypes(bxd.genotypes)

  bxd.genotypes[1:5, 1:10]
  bxd.numeric[1:5, 1:10]
  bxd.maf.genotypes[1:5, 1:10]

  require("BXDtools")
  data('bxd.phenotypes', package='BXDtools', envir=environment())
  bxd.phenotypes <- get("bxd.phenotypes", envir=environment()) 

  
  #bxd.pheno <- download.BXD.phenotypes()
  #bxd.phenotypes <- as.phenotype.matrix(bxd.genotypes, bxd.pheno)
  require("BXDtools")
  bxd.phenosomes <- only.phenosomes(bxd.phenotypes)
  scores <- do.BXD.phewas(bxd.genotypes, bxd.phenosomes, marker = "Affy_5122277", LRS = TRUE)
  plot.phewas(scores, bxd.phenosomes, TRUE, main = "BXD PheWAS at Affy_5122277 (Chr 8 - 54,959,871)")

  map <- calculate.cM.positions(bxd.genotypes)
  map[1:10,]
  
  
  require(BXDtools)
  res <- calc.qtl(bxd.genotypes, bxd.phenotypes)
  plot.qtl(bxd.genotypes, res)
}
