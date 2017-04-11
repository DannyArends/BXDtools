#
# BXD-genotypes.R
#
# copyright (c) 2017-2020 - Danny Arends, Pjotr Prins, Rob Williams, Gudrun A. Brockmann
# last modified Apr, 2017
# first written Jan, 2017
# 
# Routines to download and process BXD genotype data using GNapi
#

# Download the full BXD.geno file from genenetwork.org
download.BXD.genotypes <- function() {
  bxd.geno <- read.table("http://www.genenetwork.org/genotypes/BXD.geno", sep = "\t",
                         skip = 20, header = TRUE, row.names = 2, colClasses = "character", stringsAsFactors = FALSE)
  return(bxd.geno)
}

# Return a bxd.genotypes object without the rowheaders
only.BXD.genotypes <- function(bxd.geno) {
  BXDdatacols <- grep("BXD", colnames(bxd.geno))
  bxd.genotypes <- bxd.geno[, BXDdatacols]
  annotationcols <- seq(1:ncol(bxd.geno))[-BXDdatacols]
  attr(bxd.genotypes, "map") <- bxd.geno[, annotationcols]
  return(bxd.genotypes)
}

# Return a bxd.genotypes recoded as B = -1, H = 0, D = 1 and U = 9
recode.BXD.genotypes <- function(bxd.genotypes, coding = rbind(c("D", "1"), c("B", "-1"), c("H", "0"), c("U", "NA"))) {
  for(x in 1:nrow(coding)){
    bxd.genotypes[bxd.genotypes == coding[x, 1]] <- coding[x, 2]
  }
  bxd.numeric <- apply(bxd.genotypes, 2, function(x){
    return(as.numeric(x))
  })
  rownames(bxd.numeric) <- rownames(bxd.genotypes)
  attr(bxd.numeric, "map") <- attr(bxd.genotypes, "map")
  return(bxd.numeric)
}

# Return a bxd.genotypes recoded based on minor allele frequency (0 = Minor, 0.5 = Heterozygous, 1 = Major)
recode.BXD.to.MAF <- function(bxd.genotypes) {
  bxd.maf.genotypes <- apply(bxd.genotypes, 1, function(x) {
    tr <- table(x)[c("B", "D")]
    minorallele <- names(tr)[which.min(tr)]
    majorallele <- names(tr)[which.max(tr)]
    if(minorallele == majorallele){
      minorallele <- "B"
      majorallele <- "D"
    }
    copy.x <- x
    copy.x[copy.x == minorallele] <- 0
    copy.x[copy.x == "H"] <- 0.5
    copy.x[copy.x == "U"] <- NA
    copy.x[copy.x == majorallele] <- 1
    return(copy.x)
  })
  bxd.maf.numeric <- apply(bxd.maf.genotypes, 1, function(x) {
    return(as.numeric(x))
  })
  rownames(bxd.maf.numeric) <- rownames(bxd.genotypes)
  attr(bxd.maf.numeric, "map") <- attr(bxd.genotypes, "map")
  return(bxd.maf.numeric)
}

# Add the parental and F1 individual genotypes to the BXD genotypes
add.PandF1 <- function(bxd.genotypes) {
  if(all(unlist(apply(bxd.genotypes, 2, class)) == "numeric")) {
    # Add the parentals and F1 to a numeric matrix
    combined <- as.matrix(cbind( "B6D2F1" = rep(0, nrow(bxd.genotypes)), 
                                 "D2B6F1" = rep(0, nrow(bxd.genotypes)), 
                                 "C57BL/6J" = rep(-1, nrow(bxd.genotypes)), 
                                 "DBA/2J" = rep(1, nrow(bxd.genotypes)), 
                                 bxd.genotypes))
  } else if (all(unlist(apply(bxd.genotypes, 2, class)) == "character")) {
    # Add the parentals and F1 to a character matrix
    combined <- as.matrix(cbind( "B6D2F1" = rep("H", nrow(bxd.genotypes)), 
                                 "D2B6F1" = rep("H", nrow(bxd.genotypes)), 
                                 "C57BL/6J" = rep("B", nrow(bxd.genotypes)), 
                                 "DBA/2J" = rep("D", nrow(bxd.genotypes)), 
                                 bxd.genotypes))
  } else {
    stop("Input needs to be a numeric or character matrix")
  }
  attr(combined, "map") <- attr(bxd.genotypes, "map")
  return(combined)
}

# Get the map from a BXD genotype matrix
get.map <- function(bxd.genotypes) { return(attr(bxd.genotypes, "map")) }
