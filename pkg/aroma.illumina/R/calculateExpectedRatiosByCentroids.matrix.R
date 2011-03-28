setMethodS3("calculateExpectedRatiosByCentroids", "matrix", function(theta, muR, muT=NULL, ...) {
  # Argument 'theta':
  stopifnot(is.matrix(theta));
  dim <- dim(theta);

  # Argument 'muR':
  stopifnot(is.matrix(muR));
  stopifnot(nrow(muR) == dim[1]);
  nbrOfCentroids <- ncol(muR);
  stopifnot(nbrOfCentroids == 3);


  # Argument 'muT':
  if (!is.null(muT)) {
    stopifnot(is.matrix(muT));
    stopifnot(dim(muT) == dim(muR));
  } else {
    muT <- matrix(0, nrow=nrow(muR), ncol=nbrOfCentroids);
    muT[,2] <- 1/2;
    muT[,3] <- 1;
  }


  naValue <- as.double(NA);
  gamma <- array(naValue, dim=dim, dimnames=dimnames(theta));

  # For each sample
  for (ii in seq(length=ncol(gamma))) {
    keep <- (theta[,ii] < muT[,2]);
    keep <- keep & is.finite(keep);
    rho <- (theta[,ii] - muT[,1])  / (muT[,2]-muT[,1]);
    sR <- (muR[,2]-muR[,1]);
    Rhat <- muR[,1] + rho * sR;
    gamma[keep,ii] <- Rhat[keep];

    keep <- !keep;
    rho <- (theta[,ii] - muT[,2])  / (muT[,3]-muT[,2]);
    sR <- (muR[,3]-muR[,2]);
    Rhat <- muR[,2] + rho * sR;
    gamma[keep,ii] <- Rhat[keep];
  }

  gamma;
}) # calculateExpectedRatiosByCentroids()



############################################################################
# HISTORY:
# 2011-03-27
# o Created.
############################################################################
