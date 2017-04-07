
# Download a single BXD phenotype
download.single.BXD.phenotype <- function(datasetid, verbose = FALSE){
  info <- info_dataset(datasetid)
  if("dataset" %in% names(info) && info$dataset == "phenotype"){
    if(verbose) cat(paste0("Retrieved data for '", datasetid, "'\n"))
    return(list(info, get_pheno(datasetid)))
  }else{
    if(verbose) cat(datasetid, "not a phenotype\n")
  }
  return(NULL)
}

# Download all BXD phenotypes and descriptions in a list (of lists)
download.BXD.phenotypes <- function(verbose = FALSE){
  require("GNapi")
  datasets <- list_datasets("BXD")
  results <- vector("list", nrow(datasets))
  for(n in 1:nrow(datasets)){
    results[[n]] <- tryCatch(
      download.single.BXD.phenotype(datasets[n, "id"])
      , error = function(e) {
        if(verbose){
          cat(paste0(n, " ", datasets[n, "id"], " - '", datasets[n, "name"], "' has some error\n"))
        }
      }
    )
  }
  return(results)
}

# Convert downloaded BXD data to a matrix in R
as.phenotype.matrix <- function(bxd.genotypes, bxd.pheno){
  hasData <- which(lapply(bxd.pheno, length) == 2)
  gnIDs <- paste0("GN_", unlist(lapply(lapply(bxd.pheno[hasData],"[[", 1),"[", "id")))
  samples <- colnames(bxd.genotypes)
  phenotypes <- matrix(NA, length(hasData), length(samples), dimnames=list(gnIDs, samples))
  phenotype.descriptions <- matrix(NA, length(hasData), 4, dimnames=list(gnIDs, c("name", "title", "description", "year")))
  for(x in hasData) {
    gnid <- paste0("GN_", bxd.pheno[[x]][[1]]$id)
    gnname <- bxd.pheno[[x]][[1]]$name
    gntitle <- bxd.pheno[[x]][[1]]$title
    gndescr <- bxd.pheno[[x]][[1]]$descr
    gnyear <- bxd.pheno[[x]][[1]]$year
    if(is.null(gnname)) gnname <- "Unknown"
    if(is.null(gntitle)) gntitle <- "Unknown"
    if(is.null(gndescr)) gndescr <- "Unknown"
    if(is.null(gnyear) || gnyear == "") gnyear <- "Unknown"
    phenotype.descriptions[gnid, ] <- c(gnname, gntitle, gndescr, gnyear)
    
    trait.ind <- bxd.pheno[[x]][[2]][,"strain"]
    trait.values <- bxd.pheno[[x]][[2]][,"value"]
    names(trait.values) <- trait.ind
    
    trait.ind.with.geno <- trait.ind[which(trait.ind %in% colnames(phenotypes))]
    trait.values.with.geno <- trait.values[trait.ind.with.geno]
    
    phenotypes[gnid, trait.ind.with.geno] <- trait.values.with.geno
  }
  attr(phenotypes, "annotation") <- phenotype.descriptions
  return(phenotypes)
}

# Add a phenotype class (extracted from the description column) to each phenotype
phenotype.add.class <- function(bxd.phenotypes) {
  phenotype.descriptions <- attr(bxd.phenotypes, "annotation")
  
  splitted <- strsplit(phenotype.descriptions[,"description"], "\\,|\\:|\\;")
  classes <- unlist(lapply(splitted,"[", 1))
  classes <- gsub(" systems", "", classes)
  classes <- gsub(" system", "", classes)
  classes <- gsub(" function", "", classes)

  if(!("class" %in% colnames(phenotype.descriptions))){
    phenotype.descriptions <- cbind(phenotype.descriptions, "class" = classes)
  }else{
    phenotype.descriptions[,"class"] <- classes
  }
  attr(bxd.phenotypes, "annotation") <- phenotype.descriptions
  return(bxd.phenotypes)
}

# Get the phenotypes from phenosomes with more that a minimum number of phenotypes
only.phenosomes <- function(bxd.phenotypes, minimum = 10, verbose = FALSE) {
  bxd.phenotypes <- phenotype.add.class(bxd.phenotypes)
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

