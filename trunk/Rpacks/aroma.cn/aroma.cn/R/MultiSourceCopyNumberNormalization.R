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
#  \item{fitUgp}{An @see "AromaUgpFile" that specifies the common set of loci used to 
#    normalize the data sets at.}
#  \item{...}{Arguments passed to @see "aroma.core::AromaTabularBinaryFile".}
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
setConstructorS3("MultiSourceCopyNumberNormalization", function(dsList=NULL, fitUgp=NULL, ...) {
  if (!is.null(dsList)) {
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

    # Arguments 'fitUgp':
    className <- "AromaUgpFile";
    if (!inherits(fitUgp, className)) {
      throw("Argument 'fitUgp' is not an ", className, ": ", class(fitUgp)[1]);
    }
  }

  extend(Object(), "MultiSourceCopyNumberNormalization",
    .dsList = dsList,
    .fitUgp = fitUgp,
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
  }

  verbose && exit(verbose);

  dsSmoothList;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getSmoothedDataFiles
#
# @title "Gets a list of smoothed arrays across data sets for a particular sample"
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
setMethodS3("getSmoothedDataFiles", "MultiSourceCopyNumberNormalization", function(this, name, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  # Argument 'name':
  name <- Arguments$getCharacter(name);



  verbose && enter(verbose, "Getting list of smoothed data files for one sample");
  verbose && cat(verbose, "Sample name: ", name);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Smooth data towards target UGP, which specifies the common set of loci
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsSmoothList <- getSmoothedDataSets(this, verbose=less(verbose, 1));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data files to be used for fitting
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dfList <- sapply(dsSmoothList, function(ds) {
    idx <- indexOf(ds, name);
    df <- NULL;
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
  nok <- sapply(dfList, is.null);
  dfList <- dfList[!nok];
  rm(nok);
  verbose && cat(verbose, "Number of arrays: ", length(dfList));

  verbose && exit(verbose);

  dfList;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getSmoothedDataFiles
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
    verbose && enter(verbose, "Identify subset of (smoothing) units for fitting the model");

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
    fitUgp = getFitAromaUgpFile(this, ...)
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
setMethodS3("fitOne", "MultiSourceCopyNumberNormalization", function(this, name, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  # Argument 'name':
  name <- Arguments$getCharacter(name);

  # Argument 'force':
  force <- Arguments$getLogical(force);



  verbose && enter(verbose, "Fitting one sample across multiple sources");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get model parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  params <- getParameters(this, verbose=less(verbose, 1));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify list of data files to fit model to
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dfList <- getSmoothedDataFiles(this, name=name, verbose=less(verbose, 1));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already fitted?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fullnames <- sapply(dfList, getFullName);
  fullnames <- unname(fullnames);

  chipTypes <- sapply(dfList, getChipType);
  chipTypes <- unname(chipTypes);

  checkSums <- sapply(dfList, getChecksum);
  checkSums <- unname(checkSums);

  key <- list(method="fitOne", class="MultiSourceCopyNumberNormalization", 
              fullnames=fullnames, chipTypes=chipTypes, checkSums=checkSums);
  dirs <- c("aroma.affymetrix", "MultiSourceCopyNumberNormalization");
  if (!force) {
    transforms <- loadCache(key=key, dirs=dirs);
    if (is.null(transforms)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(transforms);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify subset of (smoothing) units for fitting the model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  units <- getSubsetToFit(this, verbose=less(verbose, 1));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract smoothed data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting data");
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);
  # Extracting data for sample to be normalized
  M <- lapply(dfList, FUN=function(df) {
    extractMatrix(df, rows=units, column=1, drop=TRUE);
  });

  rm(units);  # Not needed anymore

  M <- as.data.frame(M);
  gc <- gc();
  verbose && cat(verbose, gc);
  verbose && str(verbose, M);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit principal curve to smoothed data (M[,1], M[,2], ..., M[,K]) 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Across-source fit of smoothed data");
  t <- system.time({
    Mn <- normalizePrincipalCurve(M);
  });
  verbose && printf(verbose, "Processing time: %.1f seconds\n", as.double(t[3]));

  # Sanity check
  if (!identical(dim(Mn), dim(M))) {
    throw("Internal error: The normalize data has a different dimension that the non-normalized data: ", paste(dim(Mn), collapse="x"), " != ", paste(dim(M), collapse="x"));
  }
  gc <- gc();
  verbose && cat(verbose, gc);
  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find the normalization function for each source (array)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Finding normalization function for each array");
  transforms <- list();
  for (kk in seq(along=dfList)) {
    df <- dfList[[kk]];
    fullname <- getFullName(df);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                              kk, fullname, length(dfList)));
    transform <- makeSmoothSplinePredict(M[,kk], Mn[,kk]);

    verbose && enter(verbose, "Saving");
##      header <- list(dataSet=dsList[kk], fullname=fullnames[kk], chipType=chipTypes[kk]);
##      fit <- list(header=header, n=length(M[[kk]]), transform=transform);
##      saveObject(fit, file=pathname);
    verbose && exit(verbose);

    transforms[[kk]] <- transform;
    rm(transform);

    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  saveCache(key=key, dirs=dirs, transforms);


  transforms;
}, protected=TRUE)  # fitOne()




###########################################################################/**
# @RdocMethod fit
#
# @title "Fits the multi-source model for all samples"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments passed to @seemethod "fitOne".}
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
#   @seemethod "fitOne".
#   @seeclass
# }
#*/########################################################################### 
setMethodS3("fit", "MultiSourceCopyNumberNormalization", function(this, ..., verbose=FALSE) {
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
  verbose && enter(verbose, "Fit normalization functions to smoothed data");

  verbose && enter(verbose, "Retrieving table of samples");
  allNames <- getAllNames(this);
  nbrOfSamples <- length(allNames);
  verbose && cat(verbose, "Number of unique samples in all sets: ", nbrOfSamples);
  verbose && str(verbose, allNames);
  verbose && exit(verbose);

  units <- NULL;
  for (kk in seq(length=nbrOfSamples)) {
    name <- allNames[kk];
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", 
                                                       kk, name, nbrOfSamples));

    transforms <- fitOne(this, name=name, ..., verbose=less(verbose, 1));
    rm(transforms);

    verbose && exit(verbose);
  } # for (kk ...)

  verbose && exit(verbose);
}, protected=TRUE)






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
setMethodS3("process", "MultiSourceCopyNumberNormalization", function(this, ..., verbose=FALSE) {
})



###########################################################################
# HISTORY:
# 2008-08-18
# o Still have to write process().
# o Added Rdoc comments.
# 2008-07-04
# o Added as.character().
# o BUG FIX: getAllNames() did return duplicated names.
# 2008-06-24
# o Created first stub from existing "manual" scripts.
# 2008-05-27
# o Created "manual" script.
###########################################################################
