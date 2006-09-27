setMethodS3("getAM", "CnChipEffectFile", function(this, other, units=NULL, ...) {
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
  am;
}) # getAM()



setMethodS3("getXAM", "CnChipEffectFile", function(this, reference, chromosome, ...) {
  # Argument 'reference':
  if (!inherits(reference, "CnChipEffectFile")) {
    throw("Argument 'reference' is not an CnChipEffectFile: ", 
                                                       class(reference));
  }

  # Argument 'chromosome':
  chromosome <- Argument$getCharacter(chromosome);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve genome information, i.e. chromosome positions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving genome information");
  cdf <- getCdf(this);
  gi <- getGenomeInformation(cdf);
  units <- getUnitIndices(gi, chromosome=chromosome);
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
  am <- getAM(this, other=reference, units=units);
  nunits <- length(cn);

  cbind(x=x, am[,c("M","A")]);
}) # getXAM()


setMethodS3("fitGlad", "CnChipEffectFile", function(this, reference, chromosome, ..., verbose=FALSE) {
  require(GLAD) || throw("Package 'GLAD' not loaded.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Fitting GLAD");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached values
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="fitGlad", sample=getName(this), reference=getName(reference), chromosome=chromosome, ...);
  fit <- loadCache(key=key);
  if (!is.null(fit))
    return(fit);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit GLAD
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (X, M, A)
  verbose && enter(verbose, "Retrieving relative copy-number estimates");
  df <- getXAM(this, reference=reference, chromosome=chromosome);
  verbose && cat(verbose, sprintf("Extracted data for %d SNPs", nrow(df)));

  # Put the data in a format recognized by GLAD
  df <- data.frame(
    LogRatio=df$M, 
    PosOrder=1:nrow(df), 
    Chromosome=chromosome, 
    PosBase=df$x
  );
  verbose && str(verbose, df);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calling glad()");
  fit <- glad(df, ..., verbose=as.logical(verbose));
  verbose && exit(verbose);

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  comment <- paste(unlist(key), collapse=";");
  saveCache(fit, key=key, comment=comment);

  fit;  
}) # fitGlad()
