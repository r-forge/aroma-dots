setMethodS3("fitGlad", "CnChipEffectFile", function(this, reference, chromosome, units=NULL, ..., force=FALSE, verbose=FALSE) {
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
  key <- list(method="fitGlad", sample=getName(this), reference=getName(reference), chromosome=chromosome, units=units, ...);
  fit <- loadCache(key=key);
  if (!is.null(fit) && !force)
    return(fit);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit GLAD
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (X, M, A)
  verbose && enter(verbose, "Retrieving relative copy-number estimates");
  df <- getXAM(this, other=reference, chromosome=chromosome, units=units, verbose=less(verbose));
  verbose && str(verbose, df);
  verbose && cat(verbose, sprintf("Extracted data for %d SNPs", nrow(df)));

  # Put the data in a format recognized by GLAD
  df <- data.frame(
    LogRatio=unname(df[,"M"]), 
    PosOrder=1:nrow(df), 
    Chromosome=rep(chromosome, nrow(df)), 
    PosBase=unname(df[,"x"])
  );
  verbose && str(verbose, df);
  df <- as.profileCGH(df);
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
