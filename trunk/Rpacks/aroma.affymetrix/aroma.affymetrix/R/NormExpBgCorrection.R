###########################################################################/**
# @RdocClass NormExpBgCorrection
#
# @title "The NormExpBgCorrection class"
#
# \description{
#  @classhierarchy
#
#  This class represents the normal+exponential "background" adjustment
#  function.
#
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#     @see "BackgroundCorrection".}
#   \item{flavor}{A @character string specifying what algorithm to use.}
#   \item{args}{A @list of additional arguments passed to the 
#     correction algorithm.}
#   \item{addJitter}{If @TRUE, Zero-mean gaussian noise is added to the
#     signals before being background corrected.}
#   \item{jitterSd}{Standard deviation of the jitter noise added.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Jitter noise}{
#   The fitting algorithm of the normal+exponentital background correction
#   model may not converge if there too many small and discrete signals.  
#   To overcome this problem, a small amount of noise may be added to the 
#   signals before fitting the model.  This is an ad hoc solution that 
#   seems to work.
#   However, adding Gaussian noise may generate non-positive signals.
# }
#
# \details{
#   By default, only PM signals are background corrected and MMs are 
#   left unchanged.
# }
#
# \author{Henrik Bengtsson. 
#         Adopted from RmaBackgroundCorrection by Ken Simpson.}
#
# \seealso{
#   Internally, @see "affy::bg.adjust" is used when \code{flavor="affy"},
#   and @see "limma::backgroundCorrect" is used when \code{flavor="limma"}.
# }
#
#*/###########################################################################
setConstructorS3("NormExpBgCorrection", function(..., flavor=c("affy", "limma"), args=NULL, addJitter=FALSE, jitterSd=0.2) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'addJitter':
  addJitter <- Arguments$getLogical(addJitter);

  # Argument 'jitterSd':
  jitterSd <- Arguments$getDouble(jitterSd);

  extend(BackgroundCorrection(..., typesToUpdate="pm"), "NormExpBgCorrection",
    .flavor = flavor,
    .addJitter = addJitter,
    .jitterSd = jitterSd,
    .args = args
  );
})


setMethodS3("getParameters", "NormExpBgCorrection", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  pmOnly <- (this$.typesToUpdate == "pm");
  
  # Get parameters of this class
  params2 <- list(
    flavor = this$.flavor,
    addJitter = this$.addJitter,
    jitterSd = this$.jitterSd,
    pmOnly = pmOnly
  );

  # Algorithm parameters
  params$args <- this$.args;

  # Append the two sets
  params <- c(params, params2);

  params;
}, private=TRUE)



###########################################################################/**
# @RdocMethod process
#
# @title "Performs background correction"
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
#  Returns a @double @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "NormExpBgCorrection", function(this, ..., force=FALSE, verbose=FALSE) {

  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Background correcting data set");

  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already background corrected");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get the CDF
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf);

  # Get algorithm parameters
  params <- getParameters(this);

  # Get the output path
  outputPath <- getPath(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Try to load the require package
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  flavor <- params$flavor;
  if (flavor == "affy") {
    require("affy") || throw("Package not loaded: affy");
  } else if (flavor == "limma") {
    require("limma") || throw("Package not loaded: limma");
  }

  # Generate random jitter?
  if (params$addJitter) {
    set.seed(6022007);
    jitter <- rnorm(length(y), mean=0, sd=params$jitterSd);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # apply normal+exponential model to each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- nbrOfArrays(ds);
  verbose && enter(verbose, "Adjusting ", nbrOfArrays, " arrays");
  dataFiles <- list();
  for (kk in seq(ds)) {
    verbose && enter(verbose, sprintf("Array #%d of %d", kk, nbrOfArrays));
    df <- getFile(ds, kk);
    verbose && print(verbose, df);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Generating output pathname
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    fullname <- getFullName(df);
    filename <- sprintf("%s.CEL", fullname);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

    # Already correction?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Output data file already exists: ", pathname);
    } else {
      if (params$pmOnly) {
        verbose && enter(verbose, "Retrieving PM-only CDF indices");
        key <- list(method="getPmCellIndices", class=class(cdf)[1], chipType=chipType);
        dirs <- c("aroma.affymetrix", chipType);
        idxsPM <- loadCache(key=key, dirs=dirs);
        if (is.null(idxsPM)) {
          indices <- getCellIndices(cdf, useNames=FALSE, unlist=TRUE);
          idxsPM <- indices[isPm(cdf)];
          rm(indices);
          saveCache(idxsPM, key=key, dirs=dirs);
        }
        verbose && exit(verbose);
      } else {
        idxsPM <- NULL;
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Reading data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Extracting data");
      y <- extractMatrix(df, cells=idxsPM, drop=TRUE);
      verbose && str(verbose, y);
      verbose && exit(verbose);

      if (params$addJitter) {
        y <- y + jitter;
      }

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Correct data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Applying the normal+exponential probe model");
      verbose && cat(verbose, "Flavor: ", flavor);
      if (flavor == "affy") {
        # From package 'affy' (without a namespace)
        args <- c(list(y), params$args);
        verbose && str(verbose, args);
        y <- do.call("bg.adjust", args=args);
      } else if (flavor == "limma") {
        args <- c(list(y), params$args);
        verbose && str(verbose, args);
        y <- do.call("backgroundCorrect", args=args);
      }
      verbose && exit(verbose);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Storing data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing corrected data");
    
      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");
      createFrom(df, filename=pathname, path=NULL, verbose=less(verbose));
      verbose && exit(verbose);

      # Write calibrated data to file
      verbose2 <- -as.integer(verbose)-2;
str(idxsPM)
str(y)
      updateCel(pathname, indices=idxsPM, intensities=y, verbose=verbose2);
      rm(y, verbose2);

      verbose && exit(verbose);
    } # if (!force && isFile(pathname))

    # Assert validity of the calibrated data file
    dfC <- newInstance(df, pathname);
    # CDF inheritance
    setCdf(dfC, cdf);

    rm(df, dfC);

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  outputDataSet <- getOutputDataSet(this, force=TRUE);
  
  verbose && exit(verbose);

  invisible(outputDataSet);
})

############################################################################
# HISTORY:
# 2009-04-06
# o Verified that new NormExpBgCorrection and old RmaBackgroundCorrection
#   gives identical results.
# o Added support for 'affy' and 'limma' flavor.
# o Added support to specify and pass any algorithm parameters.
# o Created from RmaBackgroundCorrection.R.
# 2007-06-30
# o Added Rdoc comments about jitter.
# 2007-05-26
# o Updated the Rdocs.
# 2007-03-21
# o Created.
############################################################################
