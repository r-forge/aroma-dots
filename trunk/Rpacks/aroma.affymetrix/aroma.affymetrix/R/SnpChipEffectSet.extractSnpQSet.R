###########################################################################/**
# @set "class=SnpChipEffectSet"
# @RdocMethod extractSnpQSet
#
# @title "Extracts chip effects as a SnpQSet"
#
# \description{
#  @get "title".
#  The chip effects are extracted to the log (base 2) scale.
#  Note: All chip effects are copied into memory.  Moreover, any modifications
#  of the extracted data is \emph{not} reflected in the data files, and vice
#  versa.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "oligo::SnpQSet-class" object.
# }
#
# \details{
#   Simple benchmarking on shadowfax: To fit crlmm() on 90 CEPH 50K chips
#   it takes ~100 minutes, i.e. 67 s/array.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractSnpQSet", "SnpChipEffectSet", function(this, ..., force=FALSE, verbose=FALSE) {
  require(oligo) || throw("Package not loaded: oligo");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (this$mergeStrands) {
    throw("Cannot not extract SnpQSet. Strands are merged.");
  }

  verbose && enter(verbose, "Extracting chip effects as a SnpQSet for ", 
                                           nbrOfArrays(this), " arrays");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="extractSnpQSet.SnpChipEffectSet", 
                   path=getPath(this), sampleNames=getSampleNames(this));
  res <- loadCache(key=key);
  if (!force && !is.null(res)) {
    # Just a temporary fix of sample names since I've got cached data on
    # file right now. /HB 2006-10-04
    if (identical(sampleNames(res), getSampleNames(this)))
      sampleNames(res) <- getSampleNames(this);
    return(res);
  }

  cdf <- getCdf(this);
  chipType <- getChipType(cdf);
  chipType <- gsub("-monocell$", "", chipType);
  cleanChipType <- gsub("cdf$", "", cleancdfname(chipType));
  nbrOfSamples <- nbrOfFiles(this);
  snpNames <- getUnitNames(cdf);
  units <- grep("^SNP", snpNames);
  snpNames <- snpNames[units];
  fileNames <- getNames(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create necessary objects for SnpQSet
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating AnnotatedDataFrame");
  phenoData <- new("AnnotatedDataFrame", 
     data=data.frame(sample=1:nbrOfSamples, row.names=fileNames),
     varMetadata=data.frame(labelDescription="arbitrary numbering",
                                                      row.names="sample")
  );
  verbose && exit(verbose);

  verbose && enter(verbose, "Creating MIAME description");
  experimentData <- new("MIAME");
  experimentData@preprocessing$filenames <- getPathnames(this);
  experimentData@preprocessing$oligoversion <- packageDescription("oligo")$Version;
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve chip-effect estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving chip-effect estimates for ", nbrOfSamples, " arrays");
  thetas <- readUnits(this, units=units, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting by strand and allele");
  thetas <- list(
    revA=lapply(thetas, FUN=function(unit) unit[[1]][[1]]), 
    revB=lapply(thetas, FUN=function(unit) unit[[2]][[1]]), 
    fwdA=lapply(thetas, FUN=function(unit) unit[[3]][[1]]), 
    fwdB=lapply(thetas, FUN=function(unit) unit[[4]][[1]])
  );
  gc();
  verbose && exit(verbose);

  verbose && enter(verbose, "Reshaping to (log2) matrices");
  for (kk in seq(along=thetas)) {
    theta <- thetas[[kk]];
    # Extract the vector of chip effects
    theta <- unlist(theta, use.names=FALSE);
    # Chip-effects are stored on the intensity scale.  Take log2.
    theta <- log(theta, base=2);
    # Make into an JxI matrix
    theta <- matrix(theta, ncol=nbrOfSamples, byrow=TRUE);
    rownames(theta) <- snpNames; 
    colnames(theta) <- fileNames; 
    thetas[[kk]] <- theta; 
    rm(theta);
  };
  gc();
  verbose && exit(verbose);
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create SnpQSet
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Instanciating SnpQSet");
  res <- new("SnpQSet", 
    senseThetaA=thetas$fwdA, 
    senseThetaB=thetas$fwdB, 
    antisenseThetaA=thetas$revA, 
    antisenseThetaB=thetas$revB, 
    phenoData = phenoData,
    experimentData = experimentData,
    annotation = cleanChipType
  );
  rm(thetas); gc();
  verbose && exit(verbose);

  # Update the sample names
  verbose && enter(verbose, "Updating the sample names");
  sampleNames(res) <- getSampleNames(this);
  verbose && exit(verbose);

  # Save to cache!
  comment <- paste(unlist(key, use.names=FALSE), collapse=";");
  saveCache(key=key, res, comment=comment);

  verbose && exit(verbose);

  res;
})

############################################################################
# HISTORY:
# 2006-10-02
# o Added extractSnpQSet() so that we can run crlmm().
############################################################################
