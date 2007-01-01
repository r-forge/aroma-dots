###########################################################################/**
# @RdocClass PlasqModel
#
# @title "The PlasqModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the PLASQ model.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of @see "ProbePreprocessing".}
#   \item{subsetToUpdate}{The probes to be updated.
#     If @NULL, all probes are updated.}
#   \item{typesToUpdate}{Types of probes to be updated.}
#   \item{targetDistribution}{A @numeric @vector.  The empirical 
#     distribution to which all arrays should be normalized to.}
#   \item{subsetToAvg}{The probes to calculate average empirical
#     distribution over.  If a single @numeric in (0,1), then this
#     fraction of all probes will be used.  
#     If @NULL, all probes are considered.}
#   \item{typesToAvg}{Types of probes to be used when calculating the 
#     average empirical distribution.  
#     If \code{"pm"} and \code{"mm"} only perfect-match and mismatch 
#     probes are used, respectively. If \code{"pmmm"} both types are used.
#   }
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \examples{\dontrun{
#   @include "../incl/PlasqModel.Rex"
# }}
#
# @author
#*/###########################################################################
setConstructorS3("PlasqModel", function(...) {
  extend(UnitModel(...), "PlasqModel");
})


setMethodS3("getRootPath", "PlasqModel", function(this, ...) {
  "plasqData";
})

setMethodS3("getProbeParameters", "PlasqModel", function(this, ...) {
  filename <- "probeParamaters.FloatMatrix";
  pathname <- filePath(getPath(this), filename);
  if (isFile(pathname)) {
    res <- FileFloatMatrix(pathname);
  } else {
    nbrOfUnits <- nbrOfUnits(getCdf(this));
    res <- FileFloatMatrix(pathname, nrow=nbrOfUnits, ncol=12, byrow=TRUE);
  }

  res;
})

setMethodS3("getChipEffects", "PlasqModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  ces <- this$.ces;
  if (!is.null(ces))
    return(ces);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create chip-effect files 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Let the parameter object know about the CDF structure, because we 
  # might use a modified version of the one in the CEL header.
  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create chip-effect file. The CEL set is empty.");
  
  verbose && enter(verbose, "Getting copy-number chip-effect set from dataset");
  ces <- CnChipEffectSet$fromDataSet(dataset=ds, path=getPath(this), verbose=less(verbose));
  setMergeStrands(ces, TRUE);
  verbose && exit(verbose);

  # Store in cache
  this$.ces <- ces;

  ces;
})



setMethodS3("getFitUnitFunction", "PlasqModel", function(this, ...) {
  fitUnit <- function(y, type, ...) {
    # This method corresponds to EMSNP() in PLASQ500K.
    nbrOfGroups <- length(y);
  
    # Ignore non-SNP units
    if (nbrOfGroups != 2 && nbrOfGroups != 4)
      return(NULL);
    
    # Build up the data matrix and probe type vector
    yMat <- typeVec <- NULL;
    for (gg in 1:nbrOfGroups) {
      group <- .subset2(y, gg);
      value <- .subset2(group, "intensities");
      yMat <- rbind(yMat, value);
      group <- .subset2(type, gg);
      value <- .subset2(group, "plasqType");
      typeVec <- c(typeVec, value);
    }
  
    # Assert correctness
    if (nrow(yMat) != length(typeVec)) {
      throw("Internal error: Number of rows in 'yMat' and number of elements in 'typeVec' does not match: ", nrow(yMat), " != ", length(typeVec));
    }
  
    # Fit PLASQ 
#    fit <- PLASQ::EMSNP(mat=yMat, ptype=typeVec, ...);
    fit <- fitPlasqUnit(yMat, typeVec, ...);
#    fit <- EMSNP(mat=yMat, ptype=typeVec, ...);
  
    # Return
    fit;
  } # fitUnit()

  fitUnit;
})


###########################################################################/**
# @RdocMethod findUnitsTodo
#
# @title "Identifies units for which the PLM has still not been fitted to"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer @vector.
# }
#
# @author
#
# \seealso{
#   Internally this methods calls the same method for the 
#   @see "ChipEffectSet" class.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("findUnitsTodo", "PlasqModel", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  ces <- getChipEffects(this, verbose=verbose);
  findUnitsTodo(ces, ..., verbose=verbose);
})


###########################################################################/**
# @RdocMethod fit
#
# @title "Fits the model"
#
# \description{
#  @get "title" to all or to a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns the indices of the units fitted, or @NULL if no units had to 
#  be fitted.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("fit", "PlasqModel", function(this, units="remaining", ..., unitsPerChunk=moreUnits*100000/length(getDataSet(this)), moreUnits=1, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the some basic information about this model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this);
  cdf <- getCdf(ds);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  doRemaining <- FALSE;
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
  } else if (identical(units, "remaining")) {
    doRemaining <- TRUE;
  } else {
    throw("Unknown mode of argument 'units': ", mode(units));
  }

  # Argument 'unitsPerChunk':
  unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range=c(1,Inf));

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting model of class ", class(this)[1], ":");

  verbose && print(verbose, this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify units to be fitted
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
    units <- 1:nbrOfUnits;
  } else if (doRemaining) {
    verbose && enter(verbose, "Identifying non-estimated units")
    units <- findUnitsTodo(this, verbose=less(verbose));
    nbrOfUnits <- length(units);
    verbose && exit(verbose);
  } else {
    # Fit only unique units
    units <- unique(units);
    nbrOfUnits <- length(units);
  }
  verbose && printf(verbose, "Getting model fit for %d units.\n", nbrOfUnits);

  # Identify which of the requested units have *not* already been estimated
  if (!doRemaining) {
    if (force) {
      verbose && printf(verbose, "All of these are forced to be fitted.\n");
    } else {
      units <- findUnitsTodo(this, units=units, verbose=less(verbose));
      nbrOfUnits <- length(units);
      verbose && printf(verbose, "Out of these, %d units need to be fitted.\n", nbrOfUnits);
    }
  }

  # Nothing more to do?
  if (nbrOfUnits == 0)
    return(NULL);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit the model in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the model-fit function
  fitUnit <- getFitUnitFunction(this);

  # Get probe-parameter file
  paf <- getProbeParameters(this);

  # Get (and create if missing) the chip-effect files (one per array)
  ces <- getChipEffects(this, verbose=less(verbose));

  idxs <- 1:nbrOfUnits;
  head <- 1:unitsPerChunk;
  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && cat(verbose, "Number units per chunk: ", unitsPerChunk);

  # Time the fitting.
  startTime <- processTime();

  timers <- list(total=0, read=0, fit=0, writePaf=0, writeCes=0, gc=0);

  count <- 1;
  while (length(idxs) > 0) {
    tTotal <- processTime();

    verbose && enter(verbose, "Fitting chunk #", count, " of ", nbrOfChunks);
    if (length(idxs) < unitsPerChunk) {
      head <- 1:length(idxs);
    }
    uu <- idxs[head];

    verbose && cat(verbose, "Units: ");
    verbose && str(verbose, units[uu]);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get the CEL intensities by units
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nbrOfArrays <- length(ds);
    verbose && enter(verbose, "Reading probe intensities from ", nbrOfArrays, " arrays");
    tRead <- processTime();
    y <- readUnits(this, units=units[uu], ..., verbose=less(verbose));
    timers$read <- timers$read + (processTime() - tRead);
    verbose && str(verbose, y[1]);
    verbose && exit(verbose);


    verbose && enter(verbose, "Reading probe types");
    types <- readPlasqTypes(cdf, units=units[uu], ..., verbose=less(verbose));
    verbose && str(verbose, types[1]);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit the model
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Fitting unit-group model");
    tFit <- processTime();
    fit <- vector("list", length(uu));
    for (kk in seq(along=uu)) {
      fit[[kk]] <- fitUnit(.subset2(y, kk), .subset2(types, kk));
    }
    timers$fit <- timers$fit + (processTime() - tFit);
    y <- NULL; # Not needed anymore (to minimize memory usage)
    verbose && str(verbose, fit[1]);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Store probe-parameter estimates
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing probe-parameter estimates");
    tWritePaf <- processTime();
    data <- lapply(fit, FUN=function(unit) {
      coeffsF <- .subset2(unit, "coeffsF");
      coeffsR <- .subset2(unit, "coeffsR");
      sigs <- .subset2(unit, "sigs");
      c(coeffsF, sigs[1], coeffsR, sigs[2]);
    });
    data <- unlist(data, use.names=FALSE);
    data <- matrix(data, ncol=12, byrow=TRUE);
    paf[units[uu],] <- data;
    timers$writePaf <- timers$writePaf + (processTime() - tWritePaf);
    verbose && exit(verbose);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Store copy-number estimates
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing copy-number estimates");
    tWriteCes <- processTime();
    data <- lapply(fit, FUN=function(unit) {
      cn <- .subset2(unit, "CN");
      sdTheta <- rep(1, ncol(cn));
      list(
        list(theta=cn[1,], sdTheta=sdTheta),
        list(theta=cn[2,], sdTheta=sdTheta)
      )
    })
    updateUnits(ces, units=units[uu], data=data, verbose=less(verbose));
    timers$writeCes <- timers$writeCes + (processTime() - tWriteCes);
    verbose && exit(verbose);

    fit <- NULL; # Not needed anymore

    # Next chunk...
    idxs <- idxs[-head];
    count <- count + 1;

    # Garbage collection
    tGc <- processTime();
    gc();
    timers$gc <- timers$gc + (processTime() - tGc);

    timers$total <- timers$total + (processTime() - tTotal);

    verbose && exit(verbose);
  } # while()

  totalTime <- processTime() - startTime;
  if (verbose) {
    nunits <- length(units);
    t <- totalTime[3];
    printf(verbose, "Total time for all units across all %d arrays: %.2fs == %.2fmin\n", nbrOfArrays, t, t/60);
    t <- totalTime[3]/nunits
    printf(verbose, "Total time per unit across all %d arrays: %.2fs/unit\n", nbrOfArrays, t);
    t <- totalTime[3]/nunits/nbrOfArrays;
    printf(verbose, "Total time per unit and array: %.3gms/unit & array\n", 1000*t);
    t <- nbrOfUnits(cdf)*totalTime[3]/nunits/nbrOfArrays;
    printf(verbose, "Total time for one array (%d units): %.2fmin = %.2fh\n", nbrOfUnits(cdf), t/60, t/3600);
    t <- nbrOfUnits(cdf)*totalTime[3]/nunits;
    printf(verbose, "Total time for complete dataset: %.2fmin = %.2fh\n", t/60, t/3600);
    # Get distribution of what is spend where
    timers$write <- timers$writePaf + timers$writeCes;
    t <- lapply(timers, FUN=function(timer) unname(timer[3]));
    t <- unlist(t);
    t <- 100 * t / t["total"];
    printf(verbose, "Fraction of time spent on different tasks: Fitting: %.1f%%, Reading: %.1f%%, Writing: %.1f%% (of which %.2f%% is for writing chip-effects), Explicit garbage collection: %.1f%%\n", t["fit"], t["read"], t["write"], 100*t["writeCes"]/t["write"], t["gc"]);
  }

  invisible(units);
})




############################################################################
# HISTORY:
# 2007-01-01
# o Now writing probe-parameter estimates and "chip effects" to file.
# o Created.
############################################################################
