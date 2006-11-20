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
#   \item{units}{Units to be extracted.  If @NULL, all units are considered.
#    Non-SNP units are to extracted.}
#   \item{transform}{A @character string specifying what type of 
#    transformation should be applied to the chip-effect estimates, which
#    are on the intensity scale.}
#   \item{naValue}{A @double value used to replace @NA values.}
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, cached values are ignored.}
#   \item{verbose}{See @see "R.utils::Verbose".}
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
setMethodS3("extractSnpQSet", "SnpChipEffectSet", function(this, units=NULL, transform=c("log", "asinh"), naValue=-5, ..., force=FALSE, verbose=FALSE) {
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
  snpNames <- getUnitNames(cdf, units=units);
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
  missingValues <- rep(NA, nbrOfSamples);
  extractGroup <- function(unit, group) {
    if (group <= length(unit)) {
      .subset2(.subset2(unit,group),1);
    } else {
      # Handles the 500K case too.
      missingValues;
    }
  }
  thetas <- list(
    revA=lapply(thetas, FUN=extractGroup, group=1),
    revB=lapply(thetas, FUN=extractGroup, group=2),
    fwdA=lapply(thetas, FUN=extractGroup, group=3),
    fwdB=lapply(thetas, FUN=extractGroup, group=4)
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
    verbose && enter(verbose, "Group ", names(thetas)[kk]);
    theta <- thetas[[kk]];
    # Extract the vector of chip effects
    theta <- unlist(theta, use.names=FALSE);
    # Chip-effects are stored on the intensity scale.  Take log2.
    theta <- transform(theta);
    # Replace NA values
    theta[is.na(theta)] <- naValue;
    # Make into an JxI matrix
    theta <- matrix(theta, ncol=nbrOfSamples, byrow=TRUE);
    rownames(theta) <- snpNames; 
    colnames(theta) <- fileNames; 
    thetas[[kk]] <- theta; 
    rm(theta);
    verbose && exit(verbose);
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
# o BUG FIX: extractSnpQSet() would only work for Mapping50K_* chip types.
#   Code updated to work with the Mapping250K_* chip types too.
# o Added so chip effects can be extracted on the C*arcsinh scale as well 
#   as the log2 scale.
# 2006-10-02
# o Added extractSnpQSet() so that we can run crlmm().
############################################################################
