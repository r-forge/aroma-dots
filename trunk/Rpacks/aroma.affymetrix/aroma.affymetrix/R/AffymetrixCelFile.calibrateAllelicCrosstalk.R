###########################################################################/**
# @set "class=AffymetrixCelFile"
# @RdocMethod calibrateAllelicCrosstalk
#
# @title "Calibrates probepair signals for allele A allele B crosstalk"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The path where to save the calibrated data files.}
#   \item{alpha}{A @numeric @vector in [0,1] of the sequence of alpha
#     parameter in the cone fit.}
#   \item{q}{A @numeric in [0,100] of the q parameter in the cone fit.}
#   \item{Q}{A @numeric in [0,100] of the q parameter in the cone fit.}
#   \item{targetAvg}{A @numeric specifying the average signal for all 
#     probes of allele A and allele B, respectively.}
#   \item{...}{Additional arguments passed to @seemethod "getData".}
#   \item{setsOfProbes}{A named @list where each item specifies the probe indices
#     for each of the four groups; forward and reverse strands on allele A
#     and allele B.}
#   \item{overwrite}{If @TRUE, already normalized arrays are overwritten,
#     unless skipped, otherwise an error is thrown.}
#   \item{skip}{If @TRUE, the array is not normalized if it already exists.}
#   \item{format}{The output file format.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the calibrated @see "AffymetrixCelFile" object.
# }
#
# \section{Hook functions}{
#   This method call hook functions \code{*.onBegin}, \code{*.onData}, 
#   \code{*.onFit}, \code{*.onUpdated}, \code{*.onEnd}, and \code{*.onExit}
#   where \code{*} is \code{calibrateAllelicCrosstalk.AffymetrixCelFile}.
#   See code for more details.
# }
#
# @author
#
# \seealso{
#   @see "calibrateAllelicCrosstalk.AffymetrixCelSet".
#   @see "sfit::cfit".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("calibrateAllelicCrosstalk", "AffymetrixCelFile", function(this, path, alpha=c(0.1, 0.075, 0.05, 0.03, 0.01), q=2, Q=98, targetAvg=NULL, ..., setsOfProbes=NULL, overwrite=FALSE, skip=!overwrite, verbose=FALSE) {
  require(sfit) || throw("Package 'sfit' not found.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For hooks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hookName <- "calibrateAllelicCrosstalk.AffymetrixCelFile";

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  path <- Arguments$getWritablePath(path);
  if (identical(getPath(this), path)) {
    throw("Cannot calibrate data file. Argument 'path' refers to the same path as the path of the data file to be calibrated: ", path);
  }
  mkdirs(path);

  if (!is.null(targetAvg)) {
    targetAvg <- Arguments$getDouble(targetAvg, range=c(0, Inf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- sprintf("%s.cel", getName(this));
  pathname <- Arguments$getWritablePathname(filename, path=path, 
                                         mustNotExist=(!overwrite && !skip));

  # Already calibrated?
  if (isFile(pathname) && skip) {
    verbose && cat(verbose, "Calibrated data file already exists: ", pathname);
    # CDF inheritance
    res <- newInstance(this, pathname);
    setCdf(res, cdf);
    return(res);
  }

  verbose && cat(verbose, "Pathname: ", getPathname(this));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the cell indices for each possible allele basepair.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(setsOfProbes)) {
    verbose && enter(verbose, "Identifying cell indices for each possible allele basepair");
    setsOfProbes <- getAlleleProbePairs(cdf, verbose=verbose);
    gc();
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reading data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading all probe intensities");
  yAll <- getData(this, fields="intensities", ...)$intensities;
  verbose && exit(verbose);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibrating
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  callHooks(sprintf("%s.onBegin", hookName), df=this, setsOfProbes=setsOfProbes, ...);

#  library(R.graphics); 
#  Device$set(2, width="150%", height="100%");
#  subplots(length(setsOfProbes)); par(mar=c(4,4,0.5,0.5)+0.1);

  # For each allele A allele B pair...
  for (name in names(setsOfProbes)) {
    verbose && enter(verbose, "Allele basepair ", name);
    basepair <- unlist(strsplit(name, split=""));
    idx <- setsOfProbes[[name]];

    verbose && enter(verbose, "Extracting ", length(idx)/2, " probe-pair intensities");
    y <- matrix(yAll[idx], ncol=2);
    verbose && exit(verbose);

    callHooks(sprintf("%s.onData", hookName), df=this, y=y, ...);
#    lab <- basepair;
#    plot(y, pch=".", xlab=lab[1], ylab=lab[2], xlim=c(-200,65535/2), ylim=c(-200,65535/2)); abline(a=0,b=1, lty=2);

    verbose && enter(verbose, "Calibrating");

    verbose && enter(verbose, "Fitting");
    fit <- fitGenotypeCone(y, alpha=alpha, q=q, Q=Q);
    verbose && exit(verbose);
    callHooks(sprintf("%s.onFit", hookName), df=this, fit=fit, ...);

    verbose && enter(verbose, "Backtransforming");
    yC <- backtransformGenotypeCone(y, fit=fit);
    rm(fit);
    verbose && exit(verbose);

    if (!is.null(targetAvg)) {
      verbose && enter(verbose, "Rescaling to target average ", targetAvg);
      yC <- normalizeAverage(yC, targetAvg=targetAvg);
      verbose && exit(verbose);
    }

    callHooks(sprintf("%s.onUpdated", hookName), df=this, y=y, yC=yC,...);
#    points(yC, pch=".", col="red");
    rm(y);

    # Update data
    yAll[idx] <- yC;
    rm(idx, yC);

    verbose && exit(verbose);

    gc();

    verbose && exit(verbose);
  } # for (name in ...)

  rm(setsOfProbes);
  gc();

  callHooks(sprintf("%s.onEnd", hookName), df=this, yAll=yAll, ...);
#  Device$print(sprintf("%s.png", getName(this)), width="300%", height="200%");
#  dev.off();

  gc();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Storing data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Storing calibrated data");

  # Create CEL file to store results, if missing
  verbose && enter(verbose, "Creating CEL file for results, if missing");
  createFrom(this, filename=pathname, path=NULL, verbose=less(verbose));
  verbose && exit(verbose);

  # Write calibrated data to file
  verbose2 <- -as.integer(verbose)-2;
  updateCel(pathname, intensities=yAll, verbose=verbose2);
  rm(yAll);
  verbose && exit(verbose);

  # Return calibrated data file object
  res <- newInstance(this, pathname);
  # CDF inheritance
  setCdf(res, cdf);

  callHooks(sprintf("%s.onExit", hookName), df=this, dfC=res, ...);


  res;
}, private=TRUE) # calibrateAllelicCrosstalk()



############################################################################
# HISTORY:
# 2007-03-29
# o Added Rdoc comments about hook functions.
# 2006-10-06
# o Made sure CDF association is inherited.
# 2006-07-21
# o Added calibrateAllelicCrosstalk().
############################################################################
