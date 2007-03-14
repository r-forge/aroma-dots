###########################################################################/**
# @RdocClass FirmaModel
#
# @title "The FirmaModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the FIRMA (Finding Isoforms using RMA) alternative
#  splicing model.
#
# }
# 
# @synopsis 
#
# \arguments{
#   \item{rmaPlm}{An @RmaPlm object.}
#   \item{summaryMethod}{}
#   \item{tags}{}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \author{Ken Simpson (ksimpson[at]wehi.edu.au).}
#
#*/###########################################################################
setConstructorS3("FirmaModel", function(rmaPlm=NULL, summaryMethod="upperQuartile", tags="*", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    asteriskTag <- "FIRMA";

    # Update default tags
    tags[tags == "*"] <- asteriskTag;

    # Split by commas
    tags <- paste(tags, collapse=",");
    tags <- unlist(strsplit(tags, split=","));
  }

  # Argument 'summaryMethod':
  if (!is.character(summaryMethod)) {
    throw("Argument 'summaryMethod' must be a string.");
  }

  extend(UnitModel(..., tags=tags), "FirmaModel",
         .plm=rmaPlm,
         summaryMethod=summaryMethod,
         "cached:.fs"=NULL
         );
  
})


setMethodS3("getPlm", "FirmaModel", function(this, ...) {
  this$.plm;
})

setMethodS3("getDataSet", "FirmaModel", function(this, ...) {
  getDataSet(this$.plm);
})

setMethodS3("getCdf", "FirmaModel", function(this, ...) {
  getCdf(this$.plm);
})

setMethodS3("getName", "FirmaModel", function(this, ...) {
  getName(this$.plm, ...);
})

setMethodS3("getTags", "FirmaModel", function(this, ...) {
  this$.tags;
})

setMethodS3("as.character", "FirmaModel", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  ds <- getDataSet(this);
  s <- c(s, sprintf("Data set: %s", getName(ds)));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(ds))));
  s <- c(s, sprintf("Input tags: %s", getTags(this$.plm, collapse=",")));
  s <- c(s, sprintf("Output tags: %s", getTags(this, collapse=",")));
  s <- c(s, sprintf("Parameters: (%s).", getParametersAsString(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("calculateWeights", "FirmaModel", function(this, ...) {
  calculateWeights(this$.plm, ...);
})

setMethodS3("clearCache", "FirmaModel", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".paFile", ".chipFiles", ".lastPlotData")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)



setMethodS3("getFileSetClass", "FirmaModel", function(static, ...) {
  FirmaSet;
}, static=TRUE, private=TRUE)


setMethodS3("getRootPath", "FirmaModel", function(this, ...) {
  "modelFirmaModel";
}, private=TRUE)



###########################################################################/**
# @RdocMethod getFirmaScores
#
# @title "Gets the set of FIRMA results for this model"
#
# \description{
#  @get "title".
#  There is one chip-effect file per array.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "FirmaSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFirmaScores", "FirmaModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  fs <- this$.fs;
  if (!is.null(fs))
    return(fs);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create files 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Let the parameter object know about the CDF structure, because we 
  # might use a modified version of the one in the CEL header.
  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create FIRMA results file. The CEL set is empty.");
  
  verbose && enter(verbose, "Getting FIRMA results set from data set");
  # Gets the Class object
  clazz <- getFileSetClass(this);
  fs <- clazz$fromDataSet(dataSet=ds, path=getPath(this), pattern=",FIRMAscores[.](c|C)(e|E)(l|L)$", verbose=less(verbose));
  verbose && exit(verbose);

  # Store in cache
  this$.fs <- fs;

  fs;
})



###########################################################################/**
# @RdocMethod getFitFunction
#
# @title "Static method to get the low-level function that fits the PLM"
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
#  Returns a @function.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################

setMethodS3("getFitFunction", "FirmaModel", function(this, ...) {
  if (this$summaryMethod=="upperQuartile") {
    fitfcn <- function(y) {
      J <- length(y);
      list(intensities=2^(1-quantile(y, probs=0.75)), stdvs=1, pixels=1);
    }
  } else if (this$summaryMethod=="median") {
    fitfcn <- function(y) {
      J <- length(y);
      1-median(y);
      list(intensities=2^(1-median(y)), stdvs=1, pixels=1);
    }
  } else if (this$summaryMethod=="max") {
    fitfcn <- function(y) {
      J <- length(y);
      list(intensities=2^(1-max(y)), stdvs=1, pixels=1);
    }
  } else {
    fitfcn <- function(y) {
      J <- length(y);
      list(intensities=1, stdvs=1, pixels=1);
    }
  }

  fitfcn;
})

setMethodS3("getFitUnitFunction", "FirmaModel", function(this, ...) {

  fitfcn <- getFitFunction(this, ...);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Function to fit all groups (exons) within a unit
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitUnit <- function(unit, ...) {
    lapply(unit, FUN=function(group) {
      if (length(group) > 0) {
        y <- .subset2(group, 1); # Get intensities
      } else {
        y <- NULL;
      }
      y <- log(y, base=2);      # convert to log scale
      fitfcn(y);
    })
  }

  fitUnit;
}, private=TRUE)



###########################################################################/**
# @RdocMethod findUnitsTodo
#
# @title "Identifies non-fitted units"
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
#  Returns an @integer @vector of unit indices.
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
setMethodS3("findUnitsTodo", "FirmaModel", function(this, ...) {
  fs <- getFirmaScores(this);
  findUnitsTodo(fs, ...);
}, private=TRUE)



###########################################################################/**
# @RdocMethod fit
#
# @title "Estimates the model parameters"
#
# \description{
#  @get "title" for all or a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be fitted.
#     If @NULL, all units are considered.
#     If \code{remaining}, only non-fitted units are considered.
#   }
#   \item{...}{Arguments passed to @seemethod "readUnits".}
#   \item{force}{If @TRUE, already fitted units are re-fitted, and
#     cached data is re-read.}
#   \item{ram}{A @double indicating if more or less units should
#     be loaded into memory at the same time.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{moreUnits}{Deprected. Use \code{ram} instead.}
# }
#
# \value{
#  Returns an @integer @vector of indices of the units fitted, 
#  or @NULL if no units was (had to be) fitted.
# }
#
# \details{
#  All estimates are stored to file.
#
#  The non-array specific parameter estimates together with standard deviation
#  estimates and convergence information are stored in one file.
#
#  The parameter estimates specific to each array, typically "chip effects", 
#  are stored in array specific files.
#
#   Data set specific estimates [L = number of probes]:
#    phi [L doubles] (probe affinities), sd(phi) [L doubles], 
#    isOutlier(phi) [L logicals]
#
#   Algorithm-specific results:
#    iter [1 integer], convergence1 [1 logical], convergence2 [1 logical]
#    dTheta [1 double]
#    sd(eps) - [1 double] estimated standard deviation of the error term
#
#   Array-specific estimates [K = nbr of arrays]:
#    theta [K doubles] (chip effects), sd(theta) [K doubles], 
#    isOutlier(theta) [K logicals]
#   
#   For each array and each unit group, we store:
#     1 theta, 1 sd(theta), 1 isOutlier(theta), i.e. (float, float, bit)
#   => For each array and each unit (with \eqn{G_j} groups), we store:
#     \eqn{G_j} theta, \eqn{G_j} sd(theta), \eqn{G_j} isOutlier(theta), 
#   i.e. \eqn{G_j}*(float, float, bit).
#   For optimal access we store all thetas first, then all sd(theta) and the
#   all isOutlier(theta).
#   To keep track of the number of groups in each unit, we have to have a
#   (unit, ngroups) map.  This can be obtained from getUnitNames() for the
#   AffymetrixCdfFile class.
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
setMethodS3("fit", "FirmaModel", function(this, units="remaining", ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the some basic information about this model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getDataSet(this$.plm);
  cdf <- getCdf(ds);
  ws <- calculateWeights(this, verbose=verbose);
  nbrOfArrays <- length(ds);

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
  verbose && printf(verbose, "Getting FIRMA results for %d units.\n", nbrOfUnits);

  # Identify which of the requested units have *not* already been estimated
  # for all arrays
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
  # Calculate results array by array.  We should be able to process a
  # whole array in memory, so don't do in chunks for now, to speed up
  # I/O.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the model-fit function
  fitUnit <- getFitUnitFunction(this);

  # Get (and create if missing) the chip-effect files (one per array)
  fs <- getFirmaScores(this, verbose=less(verbose));

  idxs <- 1:nbrOfUnits;

  # Time the fitting.
  startTime <- processTime();

  timers <- list(total=0, read=0, fit=0, writeFs=0, gc=0);

  verbose && cat(verbose, "Units: ");
  verbose && str(verbose, units);

  unitsFull <- units;
  
  for (kk in seq(ws)) {

    units <- unitsFull;
    
    tTotal <- processTime();

    ff <- getFile(fs, kk);
    wf <- getFile(ws, kk);

    # check which units need to be fitted, for each file, to avoid having
    # to re-calculate in the event of a crash.

    verbose && enter(verbose, "Array #", kk, ": ", getName(wf));
    
    if (!force) {
      units <- findUnitsTodo(ff, units=units, force=TRUE, verbose=less(verbose));
      nbrOfUnits <- length(units);
      verbose && printf(verbose, "%d of the requested units need to be fitted.\n", nbrOfUnits);
    }

    if (length(units) > 0) {
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get the weights by unit
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      tRead <- processTime();

      y <- readUnits(wf, units=units, stratifyBy="pm", ..., force=force, verbose=less(verbose));
      timers$read <- timers$read + (processTime() - tRead);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculate FIRMA scores
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      tFit <- processTime();
      fit <- lapply(y, FUN=fitUnit);
      timers$fit <- timers$fit + (processTime() - tFit);
      y <- NULL; # Not needed anymore (to minimize memory usage)
      verbose && str(verbose, fit[1]);
      verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Store results
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing FIRMA results");
      tWriteFs <- processTime();

      suppressWarnings({
        map <- getCellMap(ff, units=units, ..., verbose=less(verbose));
      })
      intensities <- unlist(lapply(fit, function(unit) {
        lapply(unit, function(group) {
          .subset2(group, "intensities");
        })
      }), use.names=FALSE)
      stdvs <- unlist(lapply(fit, function(unit) {
        lapply(unit, function(group) {
          .subset2(group, "stdvs");
        })
      }), use.names=FALSE)
      pixels <- unlist(lapply(fit, function(unit) {
        lapply(unit, function(group) {
          .subset2(group, "pixels");
        })
      }), use.names=FALSE)
      data <- data.frame(intensities=intensities, stdvs=stdvs, pixels=pixels, cell=map[,"cell"]);
    
      updateDataFlat(ff, units=units, data=data, verbose=less(verbose));
      timers$writeFs <- timers$writeFs + (processTime() - tWriteFs);
      verbose && exit(verbose);

      fit <- NULL; # Not needed anymore
      intensities <- stdvs <- pixels <- NULL;
    
    # Garbage collection
      tGc <- processTime();
      gc <- gc();
      verbose && print(verbose, gc);
      timers$gc <- timers$gc + (processTime() - tGc);
      
      timers$total <- timers$total + (processTime() - tTotal);

    } # end of if (length(units) > 0)
    else {
      verbose && exit(verbose);
    }
      
  } # end of loop over arrays


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
    printf(verbose, "Total time for complete data set: %.2fmin = %.2fh\n", t/60, t/3600);
    # Get distribution of what is spend where
    timers$write <- timers$writeFs;
    t <- lapply(timers, FUN=function(timer) unname(timer[3]));
    t <- unlist(t);
    t <- 100 * t / t["total"];
    printf(verbose, "Fraction of time spent on different tasks: Fitting: %.1f%%, Reading: %.1f%%, Writing: %.1f%% (of which %.2f%% is for writing results), Explicit garbage collection: %.1f%%\n", t["fit"], t["read"], t["write"], 100*t["writeFs"]/t["write"], t["gc"]);
  }

  invisible(units);
})


############################################################################
# HISTORY:
# 2007-02-09
# o Created (based largely on ProbeLevelModel.R).
############################################################################
