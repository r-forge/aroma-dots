###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod bgAdjustOptical
#
# @title "Applies optical background correction to a set of CEL files"
#
# \description{
#  @get "title".
#
#  Adapted from @see "gcrma::bg.adjust.optical" in the \pkg{gcrma} package.
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The location to save the adjusted data files.}
#   \item{minimum}{The minimum adjusted intensity.  Defaults to 1.}
#   \item{subsetToUpdate}{The indices of the probes to be updated.
#     If @NULL, all are updated.}
#   \item{typesToUpdate}{Types of probes to be updated.  For more details,
#     see argument \code{types} of \code{identifyCells()} for the
#     @see "AffymetrixCdfFile" class.}
#   \item{...}{Not used.}
#   \item{overwrite}{If @TRUE, already adjusted arrays are overwritten,
#     unless skipped, otherwise an error is thrown.}
#   \item{skip}{If @TRUE, the array is not normalized if it already exists.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the background adjusted @see "AffymetrixCelFile" object.
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#*/###########################################################################
setMethodS3("bgAdjustOptical", "AffymetrixCelSet", function(this, path=NULL, name="bgOptical", subsetToUpdate=NULL, typesToUpdate=NULL, minimum=1, overwrite=FALSE, skip=!overwrite, ..., verbose=FALSE) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  if (is.null(path)) {
    # Path structure: /bgOptical/<data set name>/chip_data/<chip type>/
    path <- file.path(name, getName(this), "chip_data", getChipType(cdf));
  }
  if (!is.null(path)) {
    # Verify this path (and create if missing)
    path <- Arguments$getWritablePath(path);
  }

  if (identical(getPath(this), path)) {
    throw("Cannot calibrate data file. Argument 'path' refers to the same path as the path of the data file to be calibrated: ", path);
  }
  mkdirs(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Identifying the probes to be updated");
  subsetToUpdate <- identifyCells(cdf, indices=subsetToUpdate,
                                                     types=typesToUpdate);
  verbose && exit(verbose);

  verbose && cat(verbose, "Adjusting for optical effect for ", length(subsetToUpdate), " probes");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # optical effect correction for each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  verbose && enter(verbose, "Adjusting ", nbrOfArrays(this), " arrays");
  dataFiles <- list();
  for (kk in seq(this)) {
    verbose && enter(verbose, "Array #", kk);
    df <- getFile(this, kk);
    verbose && print(verbose, df);
    dataFiles[[kk]] <- bgAdjustOptical(df, path=path, subsetToUpdate=subsetToUpdate, typesToUpdate=NULL, minimum=1, verbose=less(verbose));
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  tmp <- newInstance(this, dataFiles);
  setCdf(tmp, getCdf(this));
  return(tmp);
  
})

###########################################################################/**
# @RdocMethod calculateGsbParameters
#
# @title "Computes parameters for adjustment of specific binding"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{nbrOfPms}{The number of random PMs to use in estimation.}
#   \item{affinities}{A @numeric @vector of probe affinities.}
#   \item{path}{If an affinities vector is not specified,
#      gives the path to a file storing the affinities.}
# }
#
# \value{
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#*/###########################################################################
setMethodS3("calculateGsbParameters", "AffymetrixCelSet", function(this, nbrOfPms=25000, affinities=NULL, path=NULL, ..., verbose=FALSE) {

  verbose <- Arguments$getVerbose(verbose);

  cdf <- getCdf(this);
  
# get path to affinities  
  if (is.null(path)) {
# try to find affinities file
    paths <- getPathnames(this)[1];
    paths <- getParent(paths);
    paths <- getParent(paths);
    paths <- paste(".",
                   paths,
                   "data/",
                   sep=";", collapse=";");

    pattern <- paste(getChipType(cdf, "-affinities.apa", sep=""));
    affinityFile <- findFiles(pattern=pattern, paths=paths, firstOnly=TRUE);
    if (is.null(affinityFile))
      throw("Could not locate probe affinities file: ", pattern);
  }

  verbose && enter(verbose, "Extracting PM indices");
  pmi <- unlist(getCellIndices(cdf, stratifyBy="pm", verbose=less(verbose, 2)));  verbose && exit(verbose);
  
  narray <- length(this);

#  set.seed(1);
#  was present in original gcrma code; left in here to allow for consistency
#  check between old and new versions

  # get a random subset of PM to use in parameter estimation
  pmi.random <- sample(pmi, nbrOfPms);
  # make sure we don't just sample from a single array; avoids problems
  # if we happened to choose a low quality or otherwise aberrant array
  iarray <- sample(1:narray, nbrOfPms, replace=TRUE);

  verbose && enter(verbose, "Extracting ", nbrOfPms, " random PM intensities from dataset");
  pm.random <- readCelIntensities(getPathnames(this), indices=pmi.random);
  verbose && exit(verbose);
  
  pm.random2 <- vector("double", 25000);

  for (i in 1:nbrOfPms) {
    pm.random2[i] <- pm.random[i, iarray[i]];
  }

  # clean up
  rm(pm.random); gc();

  verbose && enter(verbose, "Extracting probe affinities and fitting linear model")

  if (is.null(affinities)) {
    aff <- readApd(affinityFile, indices=pmi.random)$affinities;
  } else {
    aff <- affinities[pmi.random];
  }
  fit1 <- lm(log2(pm.random2) ~ aff);
  verbose && exit(verbose);
  
  return(fit1$coef);

})


###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod bgAdjustGcrma
#
# @title "Applies probe sequence based background correction to a set of
# CEL files"
#
# \description{
#  @get "title".
#
#  Adapted from @see "gcrma::bg.adjust.gcrma" in the \pkg{gcrma} package.
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The path where to save the adjusted data files.}
#   \item{name}{Name of the data set containing the background corrected
#        files.}
#   \item{type}{The type of background correction.  Currently accepted types
#       are "fullmodel" (the default, uses MMs) and "affinities" (uses
#       probe sequence only).}
#   \item{indicesNegativeControl}{Locations of any negative control
#       probes (e.g., the anti-genomic controls on the human exon array).
#       If @NULL and type=="affinities", MMs are used as the negative
#       controls.}
#   \item{opticalAdjust}{If @TRUE, apply correction for optical effect,
#       as in @see "gcrma::bg.adjust.optical".}
#   \item{gsbAdjust}{Should we adjust for specific binding (defaults to
#        @TRUE)?}
#   \item{k}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
#   \item{rho}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
#   \item{stretch}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
#   \item{fast}{If @TRUE, an ad hoc transformation of the PM is performed
#       (\code{gcrma::gcrma.bg.transformation.fast}).}
#   \item{overwrite}{If @TRUE, already adjusted arrays are overwritten,
#     unless skipped, otherwise an error is thrown.}
#   \item{skip}{If @TRUE, the array is not normalized if it already exists.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the background adjusted @see "AffymetrixCelFile" object.
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#
# \seealso{
#  @see "gcrma::bg.adjust.gcrma"
#  @seeclass
# }
#*/###########################################################################
setMethodS3("bgAdjustGcrma", "AffymetrixCelSet", function(this, path=NULL, name="bgGcrma", probePath=NULL, affinities=NULL, type="fullmodel",  indicesNegativeControl=NULL, opticalAdjust=TRUE, gsbAdjust=TRUE, k=6 * fast + 0.5 * (1 - fast), rho=0.7, stretch=1.15*fast + (1-fast), fast=TRUE, ..., verbose=FALSE) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  if (is.null(path)) {
    # Path structure: /bgGcrma/<data set name>/chip_data/<chip type>/
    path <- file.path(name, getName(this), "chip_data", getChipType(cdf));
  }
  if (!is.null(path)) {
    # Verify this path (and create if missing)
    path <- Arguments$getWritablePath(path);
  }

  if (identical(getPath(this), path)) {
    throw("Cannot calibrate data file. Argument 'path' refers to the same path as the path of the data file to be calibrated: ", path);
  }
  mkdirs(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate probe affinities, if not already existing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  if (is.null(affinities)) {
    filename <- paste(getChipType(cdf), "-affinities.apa", sep="");
    pathname <- filePath(path, filename, expandLinks="any");
    verbose && enter(verbose, "Computing probe affinities");
    if (isFile(pathname)) {
      verbose && enter(verbose, "Reading saved affinities: ", pathname);
      affinities <- readApd(pathname)$affinities;
      verbose && exit(verbose);
    } else {
      affinities <- cdf$computeAffinities(paths=probePath);
      verbose && cat(verbose, "Saving affinities: ", pathname);
      writeApd(pathname, data=affinities, name="affinities");
    }
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # optical background correction
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (opticalAdjust) {
    opticalPath <- filePath(path, "optical");
    dsOptical <- bgAdjustOptical(this, path=opticalPath, typesToUpdate="pmmm", ..., verbose=verbose);
    this <- dsOptical;
  }

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # estimate specific binding (GSB, in gcrma terminology)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  if (gsbAdjust) {
    verbose && enter(verbose, "Estimating specific binding parameters");
    gsbParameters <- calculateGsbParameters(this, affinities=affinities, path=path, ..., verbose=verbose);
    verbose && exit(verbose);
  }
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # NSB correction for each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  verbose && enter(verbose, "Adjusting ", nbrOfArrays(this), " arrays");
  dataFiles <- list();
  for (kk in seq(this)) {
    verbose && enter(verbose, "Array #", kk);
    df <- getFile(this, kk);
    verbose && print(verbose, df);
    dataFiles[[kk]] <- bgAdjustGcrma(df, path=path, gsbAdjust=gsbAdjust, gsbParameters=gsbParameters, type=type, ..., verbose=less(verbose));
    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  tmp <- newInstance(this, dataFiles);
  setCdf(tmp, getCdf(this));
  return(tmp);
  
})



############################################################################
# HISTORY:
# 2006-10-04
# o Tested, debugged, docs added.
# 2006-09-28
# o Created (based on AffymetrixCelSet.NORM.R).
############################################################################
