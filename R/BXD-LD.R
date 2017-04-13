#
# BXD-LD.R
#
# copyright (c) 2017-2020 - Danny Arends, Pjotr Prins, Rob Williams, Gudrun A. Brockmann
# last modified Apr, 2017
# first written Apr, 2017
# 
# Routines to calculate linkage disequilibrium in BXD mice
#

# Function for the optimization of the log likelihood expectation (Hardy-Weinberg) of observing P major with Q major
calc.likelihood <- function(pPQ, observed, pMajor, qMajor) {
  (2 * observed[1,1] + observed[1,2] + observed[2,1]) * log(pPQ) +                        # B with a B
  (2 * observed[1,3] + observed[1,2] + observed[2,3]) * log(pMajor-pPQ) +
  (2 * observed[3,1] + observed[2,1] + observed[3,2]) * log(qMajor-pPQ) +
  (2 * observed[3,3] + observed[3,2] + observed[2,3]) * log(1-pMajor-qMajor+pPQ) +
  observed[2,2] * log(pPQ * (1-pMajor-qMajor+pPQ) + (pMajor-pPQ) * (qMajor-pPQ))
}

# M1 and M2 coded as A, H, B
calc.LD <- function(bxd.genotypes, mname1 = "rs31443144", mname2 = "rs6269442", useEM = TRUE) {
  m1 <- bxd.genotypes[mname1, ]
  m2 <- bxd.genotypes[mname2, ]

  nodata <- c(which(is.na(m1)), which(is.na(m2)))
  if(length(nodata) > 0) {
    m1 <- m1[-nodata]
    m2 <- m2[-nodata]
  }

  # Step 1) Allele frequencies
  m1.freq <- sort(calc.alleleFrequency(m1),decreasing = TRUE)
  m2.freq <- sort(calc.alleleFrequency(m2),decreasing = TRUE)

  pMajor = m1.freq[1]; pMinor = m1.freq[2]
  qMajor = m2.freq[1]; qMinor = m2.freq[2]

  # Step 2)  Theoretical maximum range od D around: pMajor * qMajor (Dmin)
  Dmin <- max(-pMajor * qMajor, -pMinor * qMinor) #  Theoretical maximum difference Dmin (when D < 0)
  pPQmin <- pMajor * qMajor + Dmin;
  
  Dmax <- min(pMajor * qMinor, qMajor * pMinor)   #  Theoretical maximum difference Dmin (when D > 0)
  pPQmax <- pMajor * qMajor + Dmax;

  # Step 3) Observed haplotypes
  observed <- matrix(0, 3, 3, dimnames=list(c(names(m1.freq)[1], "H", names(m1.freq)[2]),c(names(m2.freq)[1], "H", names(m2.freq)[2])))
  tabular <- table(m1, m2)
  for(i in rownames(tabular)){
    for(j in colnames(tabular)){
      observed[i, j] <- tabular[i, j]
    }
  }

  # Step 4) Expectation maximization (or brute force calculation)
  if(useEM) {                                # Using the EM optimizer
    solution <- optimize(calc.likelihood, observed = observed, pMajor = pMajor, qMajor = qMajor, lower=pPQmin+.Machine$double.eps, upper=pPQmax-.Machine$double.eps, maximum = TRUE)
    pPQ <- solution$maximum
  } else {                                        # Using brute force
    s <- seq(pPQmin, pPQmax, by = 0.0001)
    lldmx <- calc.likelihood(s, observed, pMajor, qMajor)
    pPQ <- s[which.max(lldmx)]
  }

  # Step 5) Calculation of estimation of D, D' and other statistics
  estD <- pPQ - pMajor * qMajor
  if (estD > 0) {
    estDp <- estD / Dmax
  } else {
    estDp <- estD / Dmin
  }
  dchi <- as.numeric((2 * sum(observed) * estD^2) / (pMajor * pMinor * qMajor* qMinor))
  dpval <- pchisq(dchi, 1, lower.tail = FALSE)
  r <- as.numeric(estD / sqrt(pMajor * pMinor * qMajor* qMinor))

  return(list("m1.freq" = m1.freq, "m2.freq" = m2.freq , "pPQ" = pPQ, "observed" = observed, "Dmax" = Dmax, "Dmin" = Dmin, "D" = estD, "D'" = estDp, "R^2" = r^2, "sampleSize" =  sum(observed), "X^2"=dchi, "P-value" = dpval))
}

calc.LD.marker <- function(bxd.genotypes, mname = "rs31443144", distance = 10) {
  mChr <- map[mname, "Chr"]
  mPos <- as.numeric(map[mname, "Mb"])
  nearby <- which(map[,"Chr"] == mChr & as.numeric(map[, "Mb"]) > (mPos - distance) &  as.numeric(map[, "Mb"]) < (mPos + distance))
  nearby <- rownames(map)[nearby]
  LDs <- matrix(NA, length(nearby), 6)
  LDs[,1] <- rep(mname, length(nearby))
  LDs[,2] <- nearby
  pos <- 1
  for(y in nearby){
    LDres <- calc.LD(bxd.genotypes, mname, y)
    LDs[pos, 3:6] <- round(c(LDres$"m1.freq"[2], LDres$"m2.freq"[2], LDres$"D'", LDres$"R^2"), 2)
    pos <- pos + 1
  }
  colnames(LDs) <- c("From", "To", "MAF1", "MAF2", "D'", "R^2")
  return(LDs)
}

calc.LDmap <- function(bxd.genotypes, distance = 10) {
  map <- attr(bxd.genotypes, "map")
  LDmap <- NULL
  for(x in 1:nrow(map)){
    LDs <- calc.LD.marker(bxd.genotypes, rownames(map)[x], distance)
    LDmap <- rbind(LDmap, LDs)
  }
  invisible(LDmap)
}
