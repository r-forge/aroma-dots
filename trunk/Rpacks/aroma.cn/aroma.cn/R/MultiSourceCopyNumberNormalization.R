###########################################################################/**
# @RdocClass MultiSourceCopyNumberNormalization
#
# @title "The MultiSourceCopyNumberNormalization class"
#
# \description{
#  @classhierarchy
#
#  A normalization method that normalizes copy-number estimates measured
#  by multiple sites and/or platforms for common samples.  It normalizes the
#  estimates toward a common scale such that for any copy-number level 
#  the mean level of the normalized data are the same.
# }
# 
# @synopsis
#
# \arguments{
#  \item{dsList}{A @list of K @see "AromaTotalCnBinarySet":s.}
#  \item{fitUgp}{An @see "aroma.core::AromaUgpFile" that specifies the 
#    common set of loci used to normalize the data sets at.}
#  \item{subsetToFit}{The subset of loci (as mapped by the \code{fitUgp}
#    object) to be used to fit the normalization functions.
#    If @NULL, loci on chromosomes 1-22 are used, but not on ChrX and ChrY.
#  }
#  \item{targetDimension}{A @numeric index specifying the data set in
#    \code{dsList} to which each platform in standardize towards.
#    If @NULL, the arbitrary scale along the fitted principal curve
#    is used.  This always starts at zero and increases.}
#  \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   The multi-source normalization method is by nature a single-sample method,
#   that is, it normalizes arrays for one sample at the time and independently
#   of all other samples/arrays.
#
#   However, the current implementation is such that it first generates
#   smoothed data for \emph{all} samples/arrays.  Then, it normalizes the
#   sample one by one.
# }
# 
# @author
#*/###########################################################################
setConstructorS3("MultiSourceCopyNumberNormalization", function(dsList=NULL, fitUgp=NULL, subsetToFit=NULL, targetDimension=1, ...) {
  if (!is.null(dsList)) {
    # Arguments 'dsList':
    if (is.list(dsList)) {
      K <- length(dsList);

      className <- "AromaTotalCnBinarySet";
      for (kk in seq(length=K)) {
        dataSet <- dsList[[kk]];
        if (!inherits(dataSet, className)) {
          throw("Argument 'dsList' contains a non-", className, " object: ", 
                                                          class(dataSet)[1]);
        }
      }
      if (length(dsList) < 2) {
        throw("Argument 'dsList' must contain more than one ", 
                                                         className, ": ", K);
      }
    } else {
      throw("Argument 'dsList' is not a list: ", class(dsList)[1]);
    }

    # Arguments 'fitUgp':
    className <- "AromaUgpFile";
    if (!inherits(fitUgp, className)) {
      throw("Argument 'fitUgp' is not an ", className, ": ", class(fitUgp)[1]);
    }

    # Argument 'subsetToFit':
    if (is.null(subsetToFit)) {
    } else if (is.character(subsetToFit)) {
      throw("Yet not implemented: Argument 'subsetToFit' is of type character.");
    } else {
      subsetToFit <- Arguments$getIndices(subsetToFit, 
                                          range=c(1, nbrOfUnits(fitUgp)));
    }

    # Argument 'targetDimension'
    targetDimension <- Arguments$getIndex(targetDimension, range=c(1, K));
  }


  extend(Object(), "MultiSourceCopyNumberNormalization",
    .dsList = dsList,
    .fitUgp = fitUgp,
    .subsetToFit = subsetToFit,
    .targetDimension = targetDimension,
    .dsSmoothList = NULL
  )
})


setMethodS3("as.character", "MultiSourceCopyNumberNormalization", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);

  # Data sets:
  dsList <- getInputDataSets(this);
  s <- c(s, sprintf("Data sets (%d):", length(dsList)));
  for (kk in seq(along=dsList)) {
    ds <- dsList[[kk]];
    s <- c(s, as.character(ds));
  }

  # All common array names:
  names <- getAllNames(this);
  n <- length(names);
  s <- c(s, sprintf("Number of common array names: %d", n));
  if (n >= 5)
    names <- c(names[1:2], "...", names[n]);
  names <- paste(names, collapse=", ");
  s <- c(s, sprintf("Names: %s", names)); 

  # Parameters:
  paramStr <- getParametersAsString(this);
  paramStr <- paste(paramStr, collapse="; ");
  s <- c(s, sprintf("Parameters: %s", paramStr)); 

  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


###########################################################################/**
# @RdocMethod getInputDataSets
#
# @title "Gets the list of data sets to be normalized"
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
#  Returns a @list.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/########################################################################### 
setMethodS3("getInputDataSets", "MultiSourceCopyNumberNormalization", function(this, ...) {
  this$.dsList;
})



setMethodS3("nbrOfDataSets", "MultiSourceCopyNumberNormalization", function(this, ...) {
  length(getInputDataSets(this));
});



setMethodS3("getOutputPaths", "MultiSourceCopyNumberNormalization", function(this, ...) {
  dsList <- getInputDataSets(this);
  paths <- lapply(dsList, FUN=function(ds) {
    tag <- "mscn";
    path <- getPath(ds);
    path <- getParent(path, 2);
    fullname <- getFullName(ds);
    fullname <- sprintf("%s,%s", fullname, tag);
    chipType <- getChipType(ds);
    file.path(path, fullname, chipType);
  });
  paths <- unlist(paths, use.names=FALSE);
  paths;
});


setMethodS3("getOutputDataSets", "MultiSourceCopyNumberNormalization", function(this, ..., force=FALSE) {
  dsList <- getInputDataSets(this);
  paths <- getOutputPaths(this);
  dsOutList <- list();
  for (kk in seq(along=dsList)) {
    path <- paths[[kk]];
    if (isDirectory(path)) {
      ds <- dsList[[kk]];
      dsOut <- fromFiles(ds, path=path, ...);
    } else {
      dsOut <- NA;
    }
    dsOutList[[kk]] <- dsOut;
  }

  dsOutList;
})



###########################################################################/**
# @RdocMethod getFitAromaUgpFile
#
# @title "Gets the UGP file specifying the common set of loci to normalize at"
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
#  Returns a @see "aroma.core::AromaUgpFile".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/########################################################################### 
setMethodS3("getFitAromaUgpFile", "MultiSourceCopyNumberNormalization", function(this, ...) {
  this$.fitUgp;
}, protected=TRUE)




###########################################################################/**
# @RdocMethod getAllNames
#
# @title "Gets the names of all unique samples across all sources"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Passed to \code{getNames(...)} of each data set.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/########################################################################### 
setMethodS3("getAllNames", "MultiSourceCopyNumberNormalization", function(this, ...) {
  # Identify all array names across all sources
  dsList <- getInputDataSets(this);
  allNames <- lapply(dsList, getNames, ...);
  allNames <- unlist(allNames, use.names=FALSE);
  allNames <- unique(allNames);
  allNames <- sort(allNames);
  allNames;
})



###########################################################################/**
# @RdocMethod extractTupleOfDataFiles
#
# @title "Gets a list of data files for a particular name across several data sets"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{A @character string specifying the sample name of interest.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @list of K @see "AromaTotalCnBinarySet":s.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/########################################################################### 
setMethodS3("extractTupleOfDataFiles", "MultiSourceCopyNumberNormalization", function(this, dsList, name, ..., na.rm=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'dsList':
  if (is.list(dsList)) {
    className <- "AromaTotalCnBinarySet";
    for (kk in seq(along=dsList)) {
      dataSet <- dsList[[kk]];
      if (!inherits(dataSet, className)) {
        throw("Argument 'dsList' contains a non-", className, " object: ", 
                                                        class(dataSet)[1]);
      }
    }
    if (length(dsList) < 2) {
      throw("Argument 'dsList' must contain more than one ", className, 
                                                     ": ", length(dsList));
    }
  } else {
    throw("Argument 'dsList' is not a list: ", class(dsList)[1]);
  }

  # Argument 'name':
  name <- Arguments$getCharacter(name);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 



  verbose && enter(verbose, "Getting list tuple of data files for one sample");
  verbose && cat(verbose, "Sample name: ", name);

  dfList <- lapply(dsList, function(ds) {
    idx <- indexOf(ds, name);
    df <- NA;
    if (!is.na(idx)) {
      if (length(idx) > 1) {
        throw("Multiple occurances identified for this sample: ", 
                           getName(ds), " => ", paste(idx, collapse=", "));
      }
      df <- getFile(ds, idx);
    }
    df;
  });

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Filter out missing data files in order to identify the set of files
  # to fit the model on
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (na.rm) {
    keep <- sapply(dfList, FUN=function(df) !identical(df, NA));
    dfList <- dfList[keep];
  }

  verbose && cat(verbose, "Number of arrays: ", length(dfList));

  verbose && exit(verbose);

  dfList;
}, protected=TRUE)





###########################################################################/**
# @RdocMethod getSmoothedDataSets
#
# @title "Gets the data sets smoothed toward the UGP file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @list of K @see "AromaTotalCnBinarySet":s.
# }
#
# \details{
#   This method smooth each data set, each array, and each chromosome
#   toward the target (smoothing) UGP file independently of everything else.
#
#   The resulting data sets are stored in a separate location where they
#   will be located automatically in future sessions.
# }
#
# @author
#
# \seealso{
#   @seemethod "getFitAromaUgpFile".
#   @seeclass
# }
#*/########################################################################### 
setMethodS3("getSmoothedDataSets", "MultiSourceCopyNumberNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  dsSmoothList <- this$.dsSmoothList;
  if (is.null(dsSmoothList)) {
    verbose && enter(verbose, "Smoothing all data sets to the same set of loci");
    dsList <- getInputDataSets(this);
    verbose && cat(verbose, "Number of data sets: ", length(dsList));

    targetUgp <- getFitAromaUgpFile(this);
    verbose && print(verbose, targetUgp);

    sd <- 50e3;
    verbose && printf(verbose, "Kernel: %s\n", "Gaussian");
    verbose && printf(verbose, "Bandwidth (sd): %.2f\n", sd);

    dsSmoothList <- list();
    for (kk in seq(along=dsList)) {
      ds <- dsList[[kk]];
      verbose && enter(verbose, sprintf("Data set %d ('%s') of %d",
                                         kk, getFullName(ds), length(dsList)));
      sm <- TotalCnSmoothing(ds, targetUgp=targetUgp, sd=sd);
      verbose && print(verbose, sm);
      dsSmoothList[[kk]] <- process(sm, verbose=less(verbose, 1));
      verbose && exit(verbose);
    }
    names(dsSmoothList) <- names(dsList);

    # Cache in memory
    this$.dsSmoothList <- dsSmoothList;

    verbose && exit(verbose);
  }


  dsSmoothList;
}, protected=TRUE)







###########################################################################/**
# @RdocMethod getSubsetToFit
#
# @title "Gets subset of (smoothing) units for fitting the model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @integer @vector of unit indices.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/########################################################################### 
setMethodS3("getSubsetToFit", "MultiSourceCopyNumberNormalization", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  units <- this$.subsetToFit;
  if (is.null(units)) {
    verbose && enter(verbose, "Identify subset of (smoothed) units for fitting the model");

    ugp <- getFitAromaUgpFile(this);
    verbose && print(verbose, ugp);
  
    verbose && enter(verbose, "Querying UGP for units on chromosomes of interest");
    chromosomes <- 1:22;
    verbose && cat(verbose, "Chromosomes to fit: ", 
                                               seqToHumanReadable(chromosomes));
    units <- sapply(chromosomes, FUN=function(cc) {
      getUnitsOnChromosome(ugp, cc);
    });
    units <- unlist(units, use.names=FALSE);
    units <- unique(units);
    units <- sort(units);
    verbose && str(verbose, units);
    verbose && exit(verbose);

    this$.subsetToFit <- units;

    verbose && exit(verbose);
  }


  units;
}, protected=TRUE)



setMethodS3("getParameters", "MultiSourceCopyNumberNormalization", function(this, ...) {
  params <- list(
    subsetToFit = getSubsetToFit(this, ...),
    fitUgp = getFitAromaUgpFile(this, ...),
    targetDimension = this$.targetDimension
  );

  params;
})


setMethodS3("getParametersAsString", "MultiSourceCopyNumberNormalization", function(this, ...) {
  params <- getParameters(this, expand=FALSE);
  params <- trim(capture.output(str(params)))[-1];
  params <- gsub("^[$][ ]*", "", params);
  params <- gsub(" [ ]*", " ", params);
  params <- gsub("[ ]*:", ":", params);
  params;
}, private=TRUE) 



###########################################################################/**
# @RdocMethod fitOne
#
# @title "Fits the multi-source model for one sample"
#
# \description{
#  @get "title".
#  The model is fitted on the subset of units returned 
#  by @seemethod "getSubsetToFit".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{A @character string specifying the sample name of interest.}
#   \item{...}{Not used.}
#   \item{force}{If @FALSE, cached model fits are returned, otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @list of transforms.
# }
#
# @author
#
# \seealso{
#   This method is called internally by @seemethod "fit".
#   @seeclass
# }
#*/########################################################################### 
setMethodS3("fitOne", "MultiSourceCopyNumberNormalization", function(this, dfList, ..., force=FALSE, .retData=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  # Argument 'force':
  force <- Arguments$getLogical(force);



  verbose && enter(verbose, "Fitting one sample across multiple sources");
  nbrOfArrays <- length(dfList);
  verbose && cat(verbose, "Number of arrays: ", nbrOfArrays);

  # Get name of the sample from the tuple of input arrays
  # (We do it this way so that we at some stage can process() one sample
  #  at the time without first smoothing all samples. /HB 2008-08-18)
  df <- dfList[[1]];
  name <- getName(df);
  verbose && cat(verbose, "Sample name: ", name);
  rm(dfList);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get model parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- getParameters(this, verbose=less(verbose, 1));
  verbose && str(verbose, params);
  subsetToFit <- params$subsetToFit;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify list of data files to fit model to
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Smooth data towards target UGP, which specifies the common set of loci
  dsSmooth <- getSmoothedDataSets(this, verbose=less(verbose, 1));
  dfSList <- extractTupleOfDataFiles(this, dsList=dsSmooth, name=name, 
                                                 verbose=less(verbose, 1));
  rm(dsSmooth);
  verbose && str(verbose, dfSList);

  # Identify and exlude missing data sets
  keep <- sapply(dfSList, FUN=function(df) !identical(df, NA));
  keep <- whichVector(keep);
  dfSList <- dfSList[keep];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already fitted?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fullnames <- sapply(dfSList, getFullName);
  fullnames <- unname(fullnames);

  chipTypes <- sapply(dfSList, getChipType);
  chipTypes <- unname(chipTypes);

  checkSums <- sapply(dfSList, getChecksum);
  checkSums <- unname(checkSums);

  key <- list(method="fitOne", class="MultiSourceCopyNumberNormalization", 
             fullnames=fullnames, chipTypes=chipTypes, checkSums=checkSums,
             subsetToFit=subsetToFit, version="2008-10-07");
  dirs <- c("aroma.affymetrix", "MultiSourceCopyNumberNormalization");
  if (!force) {
    transforms <- loadCache(key=key, dirs=dirs);
    if (!is.null(transforms)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(transforms);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract smoothed data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting data");
  verbose && cat(verbose, "Subset of units used for fitting:");
  verbose && str(verbose, subsetToFit);
  # Extracting data for sample to be normalized
  Y <- lapply(dfSList, FUN=function(df) {
    extractMatrix(df, rows=subsetToFit, column=1, drop=TRUE);
  });

  rm(subsetToFit);  # Not needed anymore

  Y <- as.data.frame(Y);
  Y <- as.matrix(Y);
  dim <- dim(Y);
  gc <- gc();
  verbose && cat(verbose, gc);
  verbose && str(verbose, Y);
  verbose && summary(verbose, Y);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit principal curve to smoothed data (Y[,1], Y[,2], ..., Y[,K]) 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Fitting across-source normalization function");
  t <- system.time({
    fit <- fitPrincipalCurve(Y);
  });

  # Flip direction of the curve ('lambda')?
  rho <- cor(fit$lambda, Y[,1], use="complete.obs");
  flip <- (rho < 0);
  if (flip) {
    fit$lambda <- max(fit$lambda, na.rm=TRUE) - fit$lambda;
    verbose && cat(verbose, "Direction of fitted curve ('lambda') was flipped such that it increases with the signal.");
  }

  verbose && printf(verbose, "Processing time: %.1f seconds\n", 
                                                          as.double(t[3]));

  if (.retData) {
    fit$Y <- Y;
  }
  rm(Y);

  # Sanity check
  if (!identical(dim(fit$s), dim)) {
    throw("Internal error: The fitted data has a different dimension that the input data: ", paste(dim(fit$s), collapse="x"), " != ", paste(dim, collapse="x"));
  }
  gc <- gc();
  verbose && cat(verbose, gc);
  verbose && str(verbose, fit);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Standardize the channels to a target channel?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  targetChannel <- NULL;
  if (!is.null(targetChannel)) {
##     for (kk in seq(length=dim[2])) {
##       if (kk == targetChannel) {
##         targetTransform <- function(x, ...) x;
##       } else {
##         targetTransform <- makeSmoothSplinePredict(Yn[,kk], Yn[,targetChannel]);
##       }
##     } # for (kk ...)
  }

#  class(fit) <- c("MultiSourceCopyNumberNormalizationFit", class(fit));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  saveCache(key=key, dirs=dirs, fit);

  fit;
}, protected=TRUE)  # fitOne()






setMethodS3("normalizeOne", "MultiSourceCopyNumberNormalization", function(this, dfList, fit, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dfList':

  # Argument 'fit':


  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  # Argument 'force':
  force <- Arguments$getLogical(force);



  verbose && enter(verbose, "Normalize one sample across multiple sources");

  # Get name of the sample from the tuple of input arrays
  # (We do it this way so that we at some stage can process() one sample
  #  at the time without first smoothing all samples. /HB 2008-08-18)
  df <- dfList[[1]];
  name <- getName(df);
  verbose && cat(verbose, "Sample name: ", name);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get model parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- getParameters(this, verbose=less(verbose, 1));
  verbose && str(verbose, params);
  subsetToUpdate <- params$subsetToUpdate;
  targetDimension <- params$targetDimension;

  # Get (and create) the output paths
  outputPaths <- getOutputPaths(this); 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalizing array by array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing array by array");
  verbose && cat(verbose, "Units to be updated:");
  verbose && str(verbose, subsetToUpdate);

  nbrOfArrays <- length(dfList);
  dfNList <- vector("list", nbrOfArrays);
  for (kk in seq(length=nbrOfArrays)) {
    df <- dfList[[kk]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk, 
                                            getFullName(df), nbrOfArrays));

    outputPath <- outputPaths[[kk]];
    filename <- getFilename(df);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...); 
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Already normalized.");
      dfN <- newInstance(df, pathname);
    } else {
      verbose && enter(verbose, "Normalizing");

      verbose && enter(verbose, "Reading data");
      y <- extractMatrix(df, rows=subsetToUpdate, column=1, drop=TRUE);
      verbose && str(verbose, y);
      verbose && exit(verbose);
  
      verbose && enter(verbose, "Backtransforming data");
      yN <- backtransformPrincipalCurve(y, fit=fit, dimensions=kk, 
                                        targetDimension=targetDimension);
      verbose && str(verbose, yN);
      verbose && exit(verbose);

      verbose && enter(verbose, "Storing normalized data");
      verbose && cat(verbose, "Output pathname: ", pathname);

      verbose && enter(verbose, "Create output file");
      file.copy(getPathname(df), pathname);
      dfN <- newInstance(df, pathname);
      rm(pathname);
      verbose && print(verbose, dfN);
      verbose && exit(verbose);

      verbose && enter(verbose, "Writing data");
      if (is.null(subsetToUpdate)) {
        dfN[,1] <- yN;
      } else {
        dfN[subsetToUpdate,1] <- yN;
      }
      rm(yN);
      verbose && exit(verbose);

      verbose && enter(verbose, "Updating file footer");
      footer <- readFooter(dfN);
      srcFile <- df;
      footer$srcFile <- list(
        filename = getFilename(srcFile),
        filesize = getFileSize(srcFile),
        checksum = getChecksum(srcFile)
      );
      pkg <- aroma.cn;
      footer$createdBy <- list(
        class=class(this)[1],
        package = getName(pkg),
        version = getVersion(pkg)
      );
      footer$createdOn <- format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE);
      writeFooter(dfN, footer);
      verbose && exit(verbose);

      verbose && exit(verbose);

      verbose && exit(verbose);
    }

    verbose && print(verbose, dfN);
    dfNList[[kk]] <- dfN;
    rm(dfN);

    verbose && exit(verbose);
  } # for (kk ...)
  rm(subsetToUpdate);  # Not needed anymore
  gc <- gc();
  verbose && cat(verbose, gc);
  verbose && exit(verbose);

  # Return normalized arrays
  invisible(dfNList);
}, protected=TRUE)  # normalizeOne()






###########################################################################/**
# @RdocMethod process
#
# @title "Normalizes all samples"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @list of K @see "AromaTotalCnBinarySet":s.
# }
#
# @author
#
# \seealso{
#   @seemethod "fit".
#   @seeclass
# }
#*/########################################################################### 
setMethodS3("process", "MultiSourceCopyNumberNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Fit normalization functions
  # 
  # This is a multi-source (same sample across sources) whole-genome method.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Multi-source normalize all samples");
  allNames <- getAllNames(this);
  nbrOfSamples <- length(allNames);
  verbose && cat(verbose, "Number of unique samples in all sets: ", 
                                                               nbrOfSamples);
  verbose && str(verbose, allNames);

  # Get the input data sets
  dsList <- getInputDataSets(this);

  # Get (and create) the output paths
  outputPaths <- getOutputPaths(this); 


  verbose && enter(verbose, "Processing each array");
  units <- NULL;
  for (kk in seq(length=nbrOfSamples)) {
    name <- allNames[kk];
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", 
                                                    kk, name, nbrOfSamples));


    verbose && enter(verbose, "Identifying source data files");
    dfList <- extractTupleOfDataFiles(this, dsList=dsList, name=name, 
                                                   verbose=less(verbose, 1));
    verbose && print(verbose, dfList);
    verbose && exit(verbose);

   
    verbose && enter(verbose, "Check if all arrays are already normalized");
    isDone <- TRUE;
    for (jj in seq(along=dfList)) {
      df <- dfList[[jj]];
      outputPath <- outputPaths[[jj]];
      filename <- getFilename(df);
      pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...); 
      isDone <- isDone && isFile(pathname);
      if (!isDone)
        break;
    }
    verbose && cat(verbose, "Is done: ", isDone);
    verbose && exit(verbose);

    if (!force && isDone) {
      verbose && cat(verbose, "Normalized data files already exist"); 
    } else {
      verbose && enter(verbose, "Fitting model");
      fit <- fitOne(this, dfList=dfList, ..., force=force, 
                                           verbose=less(verbose, 1));
      verbose && str(verbose, fit);
      verbose && exit(verbose);
  
  
      verbose && enter(verbose, "Normalizing");
      dfNList <- normalizeOne(this, dfList=dfList, fit=fit, ..., 
                             force=force, verbose=less(verbose, 1));
      rm(fit);
      verbose && print(verbose, dfNList);

      # Sanity check
      if (length(dfNList) != length(dfList)) {
        throw("The number of normalized arrays does not match the number of source arrays: ", length(dfNList), " != ", length(dfList));
      }

      verbose && exit(verbose);
      rm(dfNList);
    }
    rm(dfList);

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  # Garbage collect
  rm(dsList);
  gc <- gc();
  verbose && print(verbose, gc);

  outputDataSets <- getOutputDataSets(this, force=TRUE, verbose=less(verbose, 1)); 

  verbose && exit(verbose);

  invisible(outputDataSets);
})


setMethodS3("pairs2", "principal.curve", function(fit, pch=19, cex=0.8, fitCol="red", fitLwd=2, fitLty=1, xlim=NULL, ylim=xlim, lower.panel=NULL, ...) {
  r <- range(c(fit$s, fit$Y), na.rm=TRUE);

  # Argument 'xlim' & 'ylim':
  if (is.null(xlim)) {
    xlim <- r;
  }
  if (is.null(ylim)) {
    ylim <- r;
  }

  hasData <- !is.null(fit$Y);

  nbrOfVars <- ncol(fit$s);
  layout(matrix(1:(nbrOfVars-1)^2, nrow=nbrOfVars-1, ncol=nbrOfVars-1, byrow=TRUE));
  par(mar=c(3,3,1,1)+0.1);
  for (rr in seq(from=1, to=nbrOfVars)) {
    if (rr == nbrOfVars)
      next;
    for (cc in seq(from=1, to=nbrOfVars)) {
      if (cc < rr) {
        plot.new();
      } else if (cc == rr) {
      } else {
        plot(NA, xlim=xlim, ylim=ylim, xlab="", ylab="");
        abline(a=0, b=1, lty=3, col="#999999", lwd=2);
        if (hasData) {
          y <- fit$Y[,c(cc,rr),drop=FALSE];
          points(y, pch=pch, cex=cex, ...);
        }

        if (fitLwd > 0) {
          y <- fit$s[,c(cc,rr),drop=FALSE];
          lines(y, col=fitCol, lwd=fitLwd, lty=fitLty);
        }
      }
    } # for (cc ...)
  } # for (rr ...)
})

###########################################################################
# HISTORY:
# 2008-10-08
# o Added argument 'targetDimension' to the constructor.
# o Now fitOne() makes sure the fitted curve has a "positive" direction.
# o Added argument 'subsetToFit' with some support, but still incomplete.
# 2008-10-07
# o Updated fitOne() and normalizeOne() to make use of the updated/new
#   fit- and backtransformPrincipalCurve() functions.
# 2008-08-18
# o Added normalizeOne() and process().
# o Added utility function extractTupleOfDataFiles().
# o Added Rdoc comments.
# 2008-07-04
# o Added as.character().
# o BUG FIX: getAllNames() did return duplicated names.
# 2008-06-24
# o Created first stub from existing "manual" scripts.
# 2008-05-27
# o Created "manual" script.
###########################################################################
