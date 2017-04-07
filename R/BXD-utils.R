
test.BXD <- function(x) {
  require("BXDtools")
  #bxd.geno <- download.BXD.genotypes()
  #bxd.genotypes <- only.BXD.genotypes(bxd.geno)
  #bxd.genotypes <- add.PandF1(bxd.genotypes)

  bxd.maf.genotypes <- recode.BXD.to.MAF(bxd.genotypes)
  bxd.numeric <- recode.BXD.genotypes(bxd.genotypes)

  bxd.genotypes[1:5, 1:10]
  bxd.numeric[1:5, 1:10]
  bxd.maf.genotypes[1:5, 1:10]

  #bxd.pheno <- download.BXD.phenotypes()
  #bxd.phenotypes <- as.phenotype.matrix(bxd.genotypes, bxd.pheno)
  bxd.phenosomes <- only.phenosomes(bxd.phenotypes)
  pvalues <- do.BXD.phewas(bxd.genotypes, bxd.phenosomes)
  plot.phewas(pvalues, bxd.phenosomes, TRUE)

  map <- calculate.cM.positions(bxd.genotypes)
  map[1:10,]
}
