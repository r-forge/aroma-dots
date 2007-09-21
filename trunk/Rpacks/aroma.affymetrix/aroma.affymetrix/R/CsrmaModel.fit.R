###########################################################################/**
# @set "class=CsrmaModel"
# @RdocMethod fit
#
# @title "Fits the CSRMA model for one chromosome across samples"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{data}{A @data.frame with columns \code{M} (log-ratio) and 
#      \code{x} (locus position).
#   }
#   \item{chromosome}{An @integer specifying the index of the chromosome to
#      be fitted.}
#   \item{...}{Additional arguments passed down to the internal fit function.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the @see "GLAD::profileCGH" object returned by @see "GLAD::glad".
# }
#
# @author
#
# \seealso{
#   Internally @see "GLAD::glad" is used.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fit", "CsrmaModel", function(this, chromosome, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'chromosome':
  chromosome <- Arguments$getIndex(chromosome);
  knownChromosomes <- getChromosomes(this);
  if (!chromosome %in% knownChromosomes) {
    throw("Argument 'chromosome' contains an unknown chromosome: ", chromosome);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting CSRMA");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- getPath(this);
  mkdirs(path);


  # Get (unit, position, chipType) map for this chromosome
  pcu <- getPositionChipTypeUnit(this, chromosome=chromosome, 
                                                  verbose=less(verbose, 10));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup data to be read
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the chip effect data sets
  cesList <- getListOfChipEffects(this, verbose=less(verbose, 10));
  



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract arguments for glad().
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- list(...);
  fitArgs <- args;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit GLAD
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up CSRMA data structure");
  nbrOfUnits <- nrow(data);
  chipTypes <- getChipTypes(this);
  data <- data.frame(
    LogRatio=data[,"M"], 
#   # Use T = M/sd (for modelling (x,T) instead. /HB 2007-02-26)
#   LogRatio=data[,"M"]/data[,"sdM"], 
    PosOrder=1:nbrOfUnits, 
    Chromosome=rep(chromosome, nbrOfUnits),
    PosBase=data[,"x"],
    # Add (chipType, units) identifiers to be able to backtrack SNP IDs etc.
    chipType=chipTypes[data[,"chipType"]],
    unit=data[,"unit"],
    # Add SD estimates
    sdTheta=data[,"sdTheta"],
    sdM=data[,"sdM"]
  );
  verbose && str(verbose, data);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calling fitCSRMA()");
  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && cat(verbose, "Chip types: ", paste(chipTypes, collapse=", "));
  verbose && cat(verbose, "Total number of units: ", nbrOfUnits);
  args <- c(list(data), fitArgs, list(verbose=as.logical(verbose)));
  rm(data, fitArgs);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # glad() writes to stdout; capture it and send it to the verbose object.
  stdout <- capture.output({
    fit <- do.call("fitCSRMA", args);
  })
  stdout <- paste(stdout, collapse="\n");
  verbose && cat(verbose, stdout);

  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;  
}, private=TRUE) # fitOne()





setMethodS3("getPositionChipTypeUnit", "CopyNumberSegmentationModel", function(this, chromosome, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'chromosome':
  chromosome <- Arguments$getIndex(chromosome);
  knownChromosomes <- getChromosomes(this);
  if (!chromosome %in% knownChromosomes) {
    throw("Argument 'chromosome' contains an unknown chromosome: ", chromosome);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (position, chipType, unit) map for this chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting (position, chipType, unit) map");

  # Get the CDFs
  cdfList <- getListOfCdfs(this, verbose=less(verbose, 10));

  # Get the genome information files  
  giList <- base::lapply(cdfList, FUN=getGenomeInformation, 
                                                  verbose=less(verbose, 10));
  verbose && print(verbose, giList);

  # Get the units on the chromosome of interest
  unitsList <- base::lapply(giList, FUN=function(gi) {
    getUnitsOnChromosome(gi, chromosome=chromosome, ...);
  });
  verbose && str(verbose, unitsList);

  # Gets (position, chipType) for these units
  posList <- vector("list", length(unitsList));
  names(posList) <- names(unitsList);
  chipTypeList <- vector("list", length(unitsList));
  names(chipTypeList) <- names(unitsList);
  for (kk in seq(along=posList)) {
    gi <- giList[[kk]];
    units <- unitsList[[kk]];
    pos <- getPositions(gi, units=units);

    # Keep only units with a position
    keep <- which(is.finite(pos));
    nbrOfUnitsBefore <- length(pos);
    nbrOfUnits <- length(keep);
    nbrOfUnitsExcl <- nbrOfUnitsBefore - nbrOfUnits;
    if (nbrOfUnitsExcl > 0) {
      pos <- pos[keep];
      units <- units[keep];
      verbose && cat(verbose, "Excluded ", nbrOfUnitsExcl, " (out of", nbrOfUnitsBefore, ") units because there is no position information available for those.");
    }
    unitsList[[kk]] <- units;
    posList[[kk]] <- pos;
    chipTypeList[[kk]] <- rep(kk, length(units));
    rm(gi, units, keep);
  }
  rm(giList);

  verbose && str(verbose, unitsList);
  verbose && str(verbose, posList);
  verbose && str(verbose, chipTypeList);

  # Unlist and order (units, position, chipType) by position
  pos <- unlist(posList, use.names=FALSE);
  rm(posList);
  o <- order(pos);
  pos <- pos[o];

  chipType <- unlist(chipTypeList, use.names=FALSE);
  rm(chipTypeList);
  chipType <- chipType[o];

  # Convert chipType into a factor
  chipTypes <- sapply(cdfList, FUN=getName);
  attr(chipType, "levels") <- chipTypes;
  class(chipType) <- "factor";

  units <- unlist(unitsList, use.names=FALSE);
  rm(unitsList);
  units <- units[o];
  rm(o);

  pcu <- data.frame(position=pos, chipType=chipType, unit=units);
  rm(units, pos, chipType);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && cat(verbose, "(position, chipType, unit) map:");
  verbose && str(verbose, pcu);

  verbose && exit(verbose);

  pcu;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2007-09-20
# o Created from GladModel.fitOne.R.
############################################################################
