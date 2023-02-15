#
# BXD-phenotypes.R
#
# copyright (c) 2017-2020 - Danny Arends, Pjotr Prins, Rob Williams, Gudrun A. Brockmann
# last modified Apr, 2017
# first written Jan, 2017
# 
# Routines to download and process BXD phenotype data using GNapi
#

# Download all BXD phenotypes and descriptions in a list (of lists)
download.BXD.phenotypes <- function(verbose = FALSE){
  tmp.annot <- tempfile()
  tmp.data <- tempfile()
  download.file("https://genenetwork.org/api/v_pre1/traits/BXDPublish.csv", tmp.annot)
  download.file("https://genenetwork.org/api/v_pre1/sample_data/BXDPublish.csv", tmp.data)
  annotation <- read.table(tmp.annot, sep = ",", header=TRUE)
  phenotypes <- read.table(tmp.data, sep = ",", header=TRUE, na.strings=c("x"), row.names=1)
  if(file.exists(tmp.annot)) file.remove(tmp.annot)
  if(file.exists(tmp.data)) file.remove(tmp.data)

  # Names should be GN_
  annotation <- annotation[!duplicated(annotation),]
  rownames(annotation) <- paste0("GN_", annotation[,"Id"])

  # Names should be GN_
  colnames(phenotypes) <- gsub("BXD_", "GN_", colnames(phenotypes))
  
  phenotypes <- phenotypes[, which(colnames(phenotypes) %in% rownames(annotation))]
  annotation <- annotation[which(rownames(annotation) %in% colnames(phenotypes)),]
  
  phenotypes <- t(phenotypes[, rownames(annotation)])
  
  phenotype.descriptions <- annotation[,c("Id", "Description", "Year")]
  colnames(phenotype.descriptions) <- c("name", "description", "year")
  rownames(phenotype.descriptions) <- rownames(annotation)
  
  phenotype.descriptions[,"description"] <- gsub("Original post publication description: ", "", phenotype.descriptions[,"description"])
  weird <- rownames(phenotype.descriptions)[grep("un-named trait", tolower(phenotype.descriptions[,"description"]))]
  phenotype.descriptions <- phenotype.descriptions[-which(rownames(phenotype.descriptions) %in% weird),]
  phenotypes <- phenotypes[-which(rownames(phenotypes) %in% weird),]
  
  attr(phenotypes, "annotation") <- phenotype.descriptions

  return(phenotypes)
}

# Convert downloaded BXD data to a matrix in R
as.phenotype.matrix <- function(bxd.genotypes, bxd.pheno){
  phenotype.descriptions <- attr(bxd.pheno, "annotation")
  bxd.pheno <- bxd.pheno[, colnames(bxd.genotypes)]
  attr(bxd.pheno, "annotation") <- phenotype.descriptions
  return(bxd.pheno)
}

# Add a phenotype class (extracted from the description column) to each phenotype
phenotype.add.class <- function(bxd.phenotypes, splitCNS = TRUE) {
  phenotype.descriptions <- attr(bxd.phenotypes, "annotation")

  splitted <- strsplit(phenotype.descriptions[,"description"], "\\,|\\:|\\;")
  classes <- tolower(unlist(lapply(splitted,"[", 1)))
  classes <- gsub(" systems", "", classes, ignore.case = TRUE)                                          # Make classes more consistent
  classes <- gsub(" system", "", classes, ignore.case = TRUE)                                           # Make classes more consistent
  classes <- gsub(" systen", "", classes, ignore.case = TRUE)                                           # Fix typos
  classes <- gsub(" sytstem", "", classes, ignore.case = TRUE)                                          # Fix typos
  classes <- gsub("enodocrine", "endocrine", classes, ignore.case = TRUE)                               # Fix typos
  classes <- gsub("endocrine", "endocrinology", classes, ignore.case = TRUE)                            # Merge endocrine into endocrinology
  classes <- gsub("reproductive", "reproduction", classes, ignore.case = TRUE)                          # Merge Reproductive into Reproduction
  classes <- gsub(" biology", "", classes, ignore.case = TRUE)                                          # Cancer biology -> Cancer
  classes <- gsub(" function", "", classes, ignore.case = TRUE)                                         # Make classes more consistent
  classes <- gsub("muscloskeletal", "musculoskeletal", classes, ignore.case = TRUE)                     # Merge Muscloskeletal into Musculoskeletal
  classes <- gsub(". metabolism", "", classes, fixed = TRUE)                                            # Fix a weird name
  classes <- gsub("central nervous", "CNS", classes, ignore.case = TRUE)                                # Rename: Central nervous to CNS
  classes <- gsub("infection disease", "infectious disease", classes, ignore.case = TRUE)               # Merge infection into infectious
  classes <- gsub(" chemistry", "", classes, ignore.case = TRUE)                                        # Make classes more consistent (blood chemistry -> blood)
  classes <- gsub("blood", "blood chemistry", classes, ignore.case = TRUE)                              # Change blood back to blood chemistry
  classes <- gsub("immune", "immune system", classes, ignore.case = TRUE)                               # Change immune back to immune system
  classes <- gsub("lung", "respiratory", classes, ignore.case = TRUE)                                   # Merge lung (2) into respiratory
  classes <- gsub("reproductive system", "reproduction", classes, ignore.case = TRUE)                   # Merge lung (2) into respiratory
  classes <- gsub("cocaine response (10 mg/kg ip)", "cocaine response", classes, fixed=TRUE)    # Merge Cocaine responses
  classes <- gsub("cocaine response (15 mg/kg ip)", "cocaine response", classes, fixed=TRUE)    # Merge Cocaine responses
  classes <- gsub("cocaine response (40 mg/kg ip)", "cocaine response", classes, fixed=TRUE)    # Merge Cocaine responses
  classes <- gsub("cocaine response (2 x 10 mg/kg ip)", "cocaine response", classes, fixed=TRUE)    # Merge Cocaine responses
  classes <- gsub("cocaine response (2 x 10 mg/kg, ip)", "cocaine response", classes, fixed=TRUE)    # Merge Cocaine responses
  classes <- gsub("cocaine response (3 x 3.2 mg/kg ip", "cocaine response", classes, fixed=TRUE)    # Merge Cocaine responses
  classes <- gsub("ethanol response (2 g/kg ip)", "ethanol response", classes, fixed=TRUE)    # Merge Ethanol responses
  classes <- gsub("ethanol response (2 g/kg)", "ethanol response", classes, fixed=TRUE)    # Merge Ethanol responses
  classes <- gsub("ethanol response (2 mg/kg ip)", "ethanol response", classes, fixed=TRUE)    # Merge Ethanol responses
  classes <- gsub("ethanol response (2.25 g/kg ip)", "ethanol response", classes, fixed=TRUE)    # Merge Ethanol responses
  classes <- gsub("ethanol response (4.0 g/kg ip)", "ethanol response", classes, fixed=TRUE)    # Merge Ethanol responses
  classes <- gsub("morphine response (50 mg/kg ip)", "morphine response", classes, fixed=TRUE)    # Merge Ethanol responses
  classes <- gsub("saline control response (10 mg/kg ip)", "saline control", classes, fixed=TRUE)    # Merge Saline responses
  classes <- gsub("saline control response (10 ml/kg ip)", "saline control", classes, fixed=TRUE)    # Merge Saline responses

  # Special handling of the CNS class since it is so big
  if(splitCNS) {
    isCNS <- which(classes == "CNS")
    classes[isCNS] <- paste0(classes[isCNS], ",", tolower(unlist(lapply(splitted[isCNS],"[", 2))))
  }

  classes <- tools::toTitleCase(classes)
  
  if(!("class" %in% colnames(phenotype.descriptions))){
    phenotype.descriptions <- cbind(phenotype.descriptions, "class" = classes)
  }else{
    phenotype.descriptions[,"class"] <- classes
  }
  attr(bxd.phenotypes, "annotation") <- phenotype.descriptions
  return(bxd.phenotypes)
}

# Get the phenotypes from phenosomes with more that a minimum number of phenotypes
only.phenosomes <- function(bxd.phenotypes, minimum = 10, splitCNS = TRUE, verbose = FALSE) {
  bxd.phenotypes <- phenotype.add.class(bxd.phenotypes, splitCNS = splitCNS)
  phenotype.descriptions <- attr(bxd.phenotypes, "annotation")
  classes <- phenotype.descriptions[, "class"]
  phenosomes <- names(which(table(classes) > minimum))
  useable.phenotypes <- which(classes %in% phenosomes)
  if(verbose) cat("Useable phenotypes:", length(useable.phenotypes), "before:", nrow(bxd.phenotypes), "\n")
  bxd.phenosomes <- bxd.phenotypes[useable.phenotypes, ]
  phenotype.descriptions <- phenotype.descriptions[useable.phenotypes, ]
  
  ordering <- sort(phenotype.descriptions[,"class"], index.return = TRUE)
  bxd.phenosomes <- bxd.phenosomes[ordering$ix, ]
  attr(bxd.phenosomes, "annotation") <- phenotype.descriptions[ordering$ix,]
  return(bxd.phenosomes)
}

# Get the map annotation from a BXD phenotype / phenosome matrix
get.annotation <- function(bxd.phenotypes) { return(attr(bxd.phenotypes, "annotation")) }

get.significant <- function(pvalues, bxd.phenosomes, alpha = 0.05, method = "BH") {
  p.adjusted <- p.adjust(pvalues, method = method)
  isSignificant <- which(p.adjusted < alpha)
  significant <- cbind(attr(bxd.phenosomes, "annotation")[isSignificant, c("name", "title", "class", "year")], Pvalue = pvalues[isSignificant], Padjusted = p.adjusted[isSignificant])
  return(data.frame(significant))
}

