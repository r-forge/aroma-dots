setConstructorS3("TotalFracBSnpData", function(data, ...) {
  colnames(data) <- c("total", "fracB");
  extend(SnpData(data, ...), "TotalFracBSnpData");
})

setMethodS3("plot", "TotalFracBSnpData", function(this, xlim=NULL, ylim=c(0,1), ...) {
  data <- this;
  NextMethod("plot", data, xlim=xlim, ylim=ylim, ...);
})


setMethodS3("asTotalFracBSnpData", "TotalFracBSnpData", function(this, ...) {
  this;
})


setMethodS3("asPolarSnpData", "TotalFracBSnpData", function(this, ...) {
  theta <- asCartesianSnpData(this);
  asPolarSnpData(theta, ...);
})

setMethodS3("asCartesianSnpData", "TotalFracBSnpData", function(this,
...) {
  totalFracB <- this;

  data <- this;
  # thetaB
  data[,2] <- totalFracB[,2]*totalFracB[,1];
  # thetaA
  data[,1] <- totalFracB[,1]-data[,2];

  CartesianSnpData(data, ...);
})


setMethodS3("callGenotypes", "TotalFracBSnpData", function(this, adjust=1.5, ...) {
  data <- this;
  nbrOfUnits <- nrow(data);

  beta <- data[,"fracB"];
  fit <- findPeaksAndValleys(beta, adjust=adjust);
  fit <- subset(fit, type == "valley");
  nbrOfGenotypeGroups <- nrow(fit)+1;
  if (nbrOfGenotypeGroups == 1) {
    warning(sprintf("PRECISION ERROR: Only one genotype group was detected for Chr%02d", chr));
    if (is.element(chr, 23:24)) {
      mu <- rep(1/2, nbrOfUnits);
      a <- 1/2;
      mu[beta < a] <- 0;
      mu[beta > a] <- 1;
    } else {
      a <- 1/3;
      b <- 2/3;
      mu[beta < a] <- 0;
      mu[beta > b] <- 1;
    }
  } else if (nbrOfGenotypeGroups == 2) {
    a <- fit$x[1];
    mu <- rep(0, nbrOfUnits);
    mu[beta > a] <- 1;
  } else if (nbrOfGenotypeGroups == 3) {
    a <- fit$x[1];
    b <- fit$x[2];
    mu <- rep(1/2, nbrOfUnits);
    mu[beta < a] <- 0;
    mu[beta > b] <- 1;
  } else {
    throw("Unexpected number of genotype groups: ", nbrOfGenotypeGroups);
  }
  mu[is.na(beta)] <- NA;

  res <- this;
  res[,2] <- mu;
  attr(res, "fit") <- fit;

  res;
})



setMethodS3("pairedBoost", "TotalFracBSnpData", function(this, dataN, ...) {
  # Coerce normal SNP signals
  tfN <- asTotalFracBSnpData(dataN);

  # Assert compatibility
  if (!identical(dim(this), dim(tfN))) {
    throw("Argument 'dataN' is of a non-compatible dimension.");
  }

  # Call normal genotypes
  tfNC <- callGenotypes(tfN, ...);

  # Estimate beta correction factors
  delta <- tfN[,2]-tfNC[,2];

  # Calibrate accordingly
  res <- this;
  res[,2] <- res[,2] - delta;
  
  res;
})


############################################################################
# HISTORY:
# 2009-03-31
# o Added pairedBoost() for TotalFracBSnpData.
# 2009-03-30
# o Created.
############################################################################
