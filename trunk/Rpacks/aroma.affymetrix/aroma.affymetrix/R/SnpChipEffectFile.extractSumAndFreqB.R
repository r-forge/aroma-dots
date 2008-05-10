setMethodS3("extractSumAndFreqB", "CnChipEffectFile", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Import thetas
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (this$combineAlleles && this$mergeStrands) {
    # theta == (theta)
    groups <- 1;
  } else if (this$combineAlleles && !this$mergeStrands) {
    # theta == (theta+, theta-)
    groups <- c(1,3);
  } else if (!this$combineAlleles && this$mergeStrands) {
    # theta == (thetaA, thetaB)
    groups <- c(1,2);
    theta <- extractThetaMatrix(this, groups=groups, ...);
    thetaSum <- rowSums(theta, na.rm=TRUE);
    rm(theta);
  } else if (!this$combineAlleles && !this$mergeStrands) {
    # theta == (thetaA+, thetaB+, thetaA-, thetaB-)
    groups <- c(1,2,3,4);
  }

  # Read data
  theta <- extractThetaMatrix(this, groups=groups, ...);
  nbrOfUnits <- nrow(theta);

  # Calculate total chip effect
  thetaSum <- rowSums(theta, na.rm=TRUE);

  # Calculate Allele B frequency?
  if (this$combineAlleles) {
    thetaRatio <- rep(NA, nbrOfUnits);
  } else {
    if (this$mergeStrands) {
      thetaRatio <- theta[,2]/thetaSum;
    } else {
      thetaRatio <- rowSums(theta[,c(2,4)], na.rm=TRUE)/thetaSum;
    }
  }
  rm(theta);

  data <- matrix(c(thetaSum, thetaRatio), nrow=nbrOfUnits, ncol=2);
  colnames(data) <- c("sum", "freqB");

  data;
})



setMethodS3("extractSumAndFreqB", "SnpChipEffectFile", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Import thetas
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (this$mergeStrands) {
    # theta == (thetaA, thetaB)
    groups <- c(1,2);
    theta <- extractThetaMatrix(this, groups=groups, ...);
    thetaSum <- rowSums(theta, na.rm=TRUE);
    rm(theta);
  } else {
    # theta == (thetaA+, thetaB+, thetaA-, thetaB-)
    groups <- c(1,2,3,4);
  }

  # Read data
  theta <- extractThetaMatrix(this, groups=groups, ...);
  nbrOfUnits <- nrow(theta);

  # Calculate total chip effect
  thetaSum <- rowSums(theta, na.rm=TRUE);

  # Calculate Allele B frequency?
  if (this$combineAlleles) {
    thetaRatio <- rep(NA, nbrOfUnits);
  } else {
    if (this$mergeStrands) {
      thetaRatio <- theta[,2]/thetaSum;
    } else {
      thetaRatio <- rowSums(theta[,c(2,4)], na.rm=TRUE)/thetaSum;
    }
  }
  rm(theta);

  data <- matrix(c(thetaSum, thetaRatio), nrow=nbrOfUnits, ncol=2);
  colnames(data) <- c("sum", "freqB");

  data;
})



############################################################################
# HISTORY:
# 2008-05-09
# o Created.
############################################################################
