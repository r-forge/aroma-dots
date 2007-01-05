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
              path=getPath(this), sampleNames=getNames(this),
              transform=transform);
  if (!force) {
    verbose && enter(verbose, "Checking cache");
    res <- loadCache(key=key);
    verbose && exit(verbose);
    if (!is.null(res))
      return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);
  chipType <- getChipType(cdf);
  chipType <- gsub("[,-]monocell$", "", chipType);
  cleanChipType <- cleancdfname(chipType, addcdf=FALSE);
  nbrOfSamples <- nbrOfFiles(this);
  units <- indexOf(cdf, pattern="^SNP");
  nbrOfUnits <- length(units);
  snpNames <- getUnitNames(cdf, units=units);
  sampleNames <- getNames(this);

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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create necessary objects for SnpQSet
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating AnnotatedDataFrame");
  phenoData <- new("AnnotatedDataFrame", 
     data=data.frame(sample=1:nbrOfSamples, row.names=sampleNames),
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

  verbose && enter(verbose, "Retrieving the target strand information for all unit groups");
  gs <- readCdfGroupStrands(getPathname(cdf), units=units, what="target");
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data by antisense and sense strands
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reordering groups by strandiness");
  verbose && cat(verbose, "Number of SNPs: ", nbrOfUnits);

  # 1) Reverse all (sense, sense, antisense, antisense) units
  value <- c("sense", "sense", "antisense", "antisense");
  idxs <- unlist(lapply(gs, FUN=identical, value), use.names=FALSE);
  thetas[idxs] <- lapply(thetas[idxs], FUN=function(groups) groups[c(3,4,1,2)]);

  # 2) Add empty (filled with NAs) groups for single-stranded SNPs
  missingValues <- rep(NA, nbrOfSamples);
  missingGroups <- list(missingValues, missingValues);

  # 2a) Add missing 'antisense' groups
  value <- c("sense", "sense");
  idxs <- unlist(lapply(gs, FUN=identical, value), use.names=FALSE);
  thetas[idxs] <- lapply(thetas[idxs], FUN=function(groups) c(missingGroups, groups));

  # 2b) Add missing 'sense' groups
  value <- c("antisense", "antisense");
  idxs <- unlist(lapply(gs, FUN=identical, value), use.names=FALSE);
  thetas[idxs] <- lapply(thetas[idxs], FUN=function(groups) c(groups, missingGroups));

  rm(gs); # Not needed anymore
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reshape and transform
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Now all units have four groups each with nbrOfSamples elements...
  verbose && enter(verbose, "Reshaping and transforming data");
  thetas <- unlist(thetas, use.names=FALSE);

  verbose && enter(verbose, "Transforming");
  # Chip-effects are stored on the intensity scale.
  thetas <- transform(thetas);
  # Replace NA values
#  thetas[is.na(thetas)] <- naValue;
  thetas[!is.finite(thetas)] <- naValue;
  verbose && exit(verbose);

  verbose && enter(verbose, "Reshaping");
  verbose && str(verbose, thetas);
  dimnames <- list(sampleNames, c("A-", "B-", "A+", "B+"), snpNames);
  dim <- c(nbrOfSamples, 4, nbrOfUnits);
  thetas <- array(thetas, dim=dim, dimnames=dimnames);
  rm(dim, dimnames);
  # Dimension order: (unit, group, sample)
  thetas <- aperm(thetas, perm=c(3,2,1));

  # Reorder SNPs in lexicographic ordering (as in oligo)
  o <- order(dimnames(thetas)[[1]]);
  thetas <- thetas[o,,];
  rm(o);
  verbose && cat(verbose, "Units reordered by lexicographic ordering:");
  verbose && str(verbose, thetas);
  verbose && exit(verbose);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create SnpQSet
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Instanciating SnpQSet");
  res <- new("SnpQSet", 
    senseThetaA = thetas[,"A+",], 
    senseThetaB = thetas[,"B+",], 
    antisenseThetaA = thetas[,"A-",], 
    antisenseThetaB = thetas[,"B-",], 
    phenoData = phenoData,
    experimentData = experimentData,
    annotation = cleanChipType
  );
  rm(thetas);
  verbose && exit(verbose);

  # Update the sample names
  verbose && enter(verbose, "Updating the sample names");
#  sampleNames(theta) <- getNames(this);
  verbose && exit(verbose);

  # Save to cache!
  comment <- paste(unlist(key, use.names=FALSE), collapse=";");
  saveCache(key=key, res, comment=comment);

  verbose && exit(verbose);

  res;
})

############################################################################
# HISTORY:
# 2006-12-12
# o Updated so that the strandiness for each SNP is correct (I think).
# 2006-10-05
# o BUG FIX: extractSnpQSet() would only work for Mapping50K_* chip types.
#   Code updated to work with the Mapping250K_* chip types too.
# o Added so chip effects can be extracted on the C*arcsinh scale as well 
#   as the log2 scale.
# 2006-10-02
# o Added extractSnpQSet() so that we can run crlmm().
############################################################################
