setMethodS3("getAM", "ChipEffectFile", function(this, other, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'other':
  if (!inherits(other, "CnChipEffectFile")) {
    throw("Argument 'other' is not an CnChipEffectFile: ", class(other));
  }

  # Argument 'units':
  cdf <- getCdf(this);
  if (is.null(units)) {
  } else {
    units <- Arguments$getIndices(units, range=c(1,nbrOfUnits(cdf)));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting (A,M)-transformed chip effects");

  nunits <- length(units);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the sample
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving sample thetas");
  theta <- unlist(this[units], use.names=FALSE);
  stopifnot(identical(length(theta), nunits));
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the other
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving other thetas");
  # Workaround for now (just in case). /HB 2006-09-26
  other$mergeStrands <- this$mergeStrands;
  other$combineAlleles <- this$combineAlleles;
  # Get the other theta estimates
  thetaRef <- unlist(other[units], use.names=FALSE);
  stopifnot(identical(length(thetaRef), nunits));
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate raw copy numbers relative to the other
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  M <- log(theta/thetaRef, base=2);
  A <- log(theta*thetaRef, base=2)/2;
  stopifnot(identical(length(M), nunits));

  am <- matrix(c(A,M), ncol=2)
  colnames(am) <- c("A", "M");
  rownames(am) <- units;

  verbose && exit(verbose);

  am;
}) # getAM()



setMethodS3("getXAM", "ChipEffectFile", function(this, other, chromosome, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'other':
  if (!inherits(other, "CnChipEffectFile")) {
    throw("Argument 'other' is not an CnChipEffectFile: ", class(other));
  }

  # Argument 'chromosome':
  chromosome <- Arguments$getCharacter(chromosome);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Getting (X,A,M)-transformed chip effects");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve genome information, i.e. chromosome positions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving genome information");
  cdf <- getCdf(this);
  gi <- getGenomeInformation(cdf);
  verbose && str(verbose, units);
  units <- getUnitIndices(gi, chromosome=chromosome, units=units, verbose=less(verbose));
  nunits <- length(units);
  if (nunits == 0)
    throw("No SNPs found on requested chromosome: ", chromosome);
  # Get the positions of all SNPs
  x <- getPositions(gi, units=units);
  stopifnot(identical(length(x), nunits));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove SNPs for which we have no position information
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  keep <- which(!is.na(x));
  nexcl <- length(x) - length(keep);
  if (nexcl > 0) {
    msg <- sprintf("Could not find position information on %d SNPs: ", nexcl);
    verbose && cat(verbose, msg);
    warning(msg);
    x <- x[keep];
    units <- units[keep];
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the relative copy-number estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  am <- getAM(this, other=other, units=units, verbose=less(verbose));
  xam <- cbind(x=x, am[,c("M","A"), drop=FALSE]);

  verbose && exit(verbose);

  xam;
}) # getXAM()


############################################################################
# HISTORY:
# 2006-09-26
# o Created.
############################################################################
