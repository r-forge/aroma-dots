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
#   Simple benchmarking on shadowfax: 
#   To extract chip effects for the 90 CEPH Mapping50K\_Xba240 arrays,
#   it takes ~90 seconds, i.e. 1s/array.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractSnpQSet", "SnpChipEffectSet", function(this, transform=c("log", "asinh"), ..., force=FALSE, verbose=FALSE) {
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

  if (is.character(transform)) {
    transform <- match.arg(transform);
  } else if (is.function(transform)) {
  } else {
    throw("Argument 'transform' must be a character of a function: ", mode(transform));
  }

  verbose && enter(verbose, "Extracting chip effects as a SnpQSet for ", 
                                           nbrOfArrays(this), " arrays");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="extractSnpQSet.SnpChipEffectSet", 
              path=getPath(this), sampleNames=getSampleNames(this),
              transform=transform);
  if (!force) {
    verbose && enter(verbose, "Checking cache");
    res <- loadCache(key=key);
    verbose && exit(verbose);
    if (!is.null(res))
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


  if (is.character(transform)) {
    if (identical(transform, "log")) {
      transform <- function(x) { log(x, base=2) }
    } else if (identical(transform, "asinh")) {
      # Constant to make arcsinh transform values to be
      # on the same scale as the log2 transformed ones.
      C <- log(2^16, base=2)/asinh(2^16);
      transform <- function(x) { C*asinh(x) }
    } else {
      throw("Unknown transform: ", transform);
    }
  }

  verbose && enter(verbose, "Reshaping to (transformed) matrices");
  for (kk in seq(along=thetas)) {
    theta <- thetas[[kk]];
    # Extract the vector of chip effects
    theta <- unlist(theta, use.names=FALSE);
    # Chip-effects are stored on the intensity scale.  Take log2.
    theta <- transform(theta);
    theta[is.na(theta)] <- -5;
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
# 2006-10-05
# o Added so chip effects can be extracted on the C*arcsinh scale as well 
#   as the log2 scale.
# 2006-10-02
# o Added extractSnpQSet() so that we can run crlmm().
############################################################################
