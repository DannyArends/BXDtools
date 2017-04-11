#
# BXD-allelefrequencies.R
#
# copyright (c) 2017-2020 - Danny Arends, Pjotr Prins, Rob Williams, Gudrun A. Brockmann
# last modified Apr, 2017
# first written Apr, 2017
# 
# Routines to calculate allele frequencies
#

calc.alleleFrequencies <- function(bxd.genotypes) {
  allelefreq <- apply(bxd.genotypes, 1, function(x){
    tbl <- table(x)
    return(c((tbl["B"] * 2 + tbl["H"]) / (2* sum(tbl)), (tbl["D"] * 2 + tbl["H"]) / (2* sum(tbl))))
  })
  allelefreq <- t(allelefreq)
  allelefreq[,1] <- round(as.numeric(allelefreq[,1]), digits = 3)
  allelefreq[,2] <- round(as.numeric(allelefreq[,2]), digits = 3)
  colnames(allelefreq) <- c("B", "D")
  return(allelefreq)
}

calc.MAF <- function(bxd.genotypes) {
  allelefreq <- apply(bxd.genotypes, 1, function(x){
    tbl <- table(x)
    if(tbl["B"] < tbl["D"])  return( c("B", (tbl["B"] * 2 + tbl["H"]) / (2* sum(tbl))) )
    if(tbl["D"] < tbl["B"])  return( c("D", (tbl["D"] * 2 + tbl["H"]) / (2* sum(tbl))) )
    if(tbl["D"] == tbl["B"]) return( c("E", 0.5) )
  })
  allelefreq <- t(allelefreq)
  allelefreq[,2] <- round(as.numeric(allelefreq[,2]), digits = 3)
  colnames(allelefreq) <- c("Allele", "MAF")
  return(allelefreq)
}
