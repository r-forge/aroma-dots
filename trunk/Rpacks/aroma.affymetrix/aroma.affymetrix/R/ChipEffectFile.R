###########################################################################/**
# @RdocClass ChipEffectFile
#
# @title "The ChipEffectFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in the probe-level models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ParameterCelFile".}
#   \item{model}{The specific type of model, e.g. \code{"pm"}.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#   An object of this class is typically obtained through the
#   \code{getChipEffects()} method for the @see "ProbeLevelModel" class.
#   An object of this class is typically part of a @see "ChipEffectSet".
# }
#*/###########################################################################
setConstructorS3("ChipEffectFile", function(..., model=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  this <- extend(ParameterCelFile(...), "ChipEffectFile",
    "cached:.firstCells" = NULL,
    model = model
  )

  setEncodeFunction(this, function(groupData, ...) {
    theta <- .subset2(groupData, "theta");
    stdvs <- .subset2(groupData, "sdTheta");
    outliers <- .subset2(groupData, "thetaOutliers");
    pixels <- NULL;
    if (!is.null(outliers))
      pixels <- -as.integer(outliers);
  
    res <- list();
    if (!is.null(theta))
      res$intensities <- theta;
    if (!is.null(stdvs))
      res$stdvs <- stdvs;
    if (!is.null(pixels))
      res$pixels <- pixels;
  
    res;
  })
  
  setDecodeFunction(this, function(groupData, ...) {
    res <- list();
    if (!is.null(groupData$intensities))
      res$theta <- groupData$intensities;
    if (!is.null(groupData$stdvs))
      res$sdTheta <- groupData$stdvs;
    if (!is.null(groupData$pixels))
      res$thetaOutliers <- as.logical(-groupData$pixels);
    res;
  })
 
  this;
})


setMethodS3("createParamCdf", "ChipEffectFile", function(static, sourceCdf, ..., verbose=FALSE) {
  # Argument 'verbose': 
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Creating CDF for chip effects");
  verbose && cat(verbose, "Source chip type: ", getChipType(sourceCdf));
  verbose && cat(verbose, "Source CDF: ", getPathname(sourceCdf));

  # CDF for chip effects
  chipType <- sprintf("%s-monocell", getChipType(sourceCdf));
  verbose && cat(verbose, "Chip type for chip effects: ", getChipType(sourceCdf));

  # Search for CDF
  cdf <- findCdf(chipType);
  if (is.null(cdf)) {
    verbose && cat(verbose, "Pathname: Not found!");
    verbose && cat(verbose, "Will create from the CDF of the dataset. NOTE: This will take several minutes or more!");
    verbose && enter(verbose, "Creating CDF");
    cdf <- createMonoCell(sourceCdf, verbose=less(verbose));
    verbose && exit(verbose);
  } else {
    verbose && cat(verbose, "Pathname: ", cdf);
    cdf <- AffymetrixCdfFile$fromFile(cdf);
  }
  verbose && exit(verbose);

  cdf;
}, static=TRUE)


setMethodS3("fromDataFile", "ChipEffectFile", function(static, df, filename=sprintf("%s-chipEffects.cel", getName(df)), path, name=getName(df), cdf=NULL, ..., verbose=FALSE) {
  # Argument 'df':
  if (!inherits(df, "AffymetrixCelFile")) {
    throw("Argument 'df' is not an AffymetrixCelFile: ", class(df)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  pathname <- Arguments$getWritablePathname(filename, path=path);
  if (!isFile(pathname)) {
    verbose && enter(verbose, "Creating chip-effect file");
    verbose && cat(verbose, "Pathname: ", pathname);

    # Get CDF for chip effects
    if (is.null(cdf))
      cdf <- createParamCdf(static, getCdf(df), verbose=less(verbose));
  
    # Get CDF header
    cdfHeader <- getHeader(cdf);

    celHeader0 <- readCelHeader(getPathname(df));
  
    # Build a valid CEL header
    celHeader <- cdfHeaderToCelHeader(cdfHeader, sampleName=name);

    # Add some extra information about what the CEL file is for
    params <- c(Descripion="This CEL file contains chip-effect estimates from the aroma.affymetrix package.");
    parameters <- gsub(" ", "_", params);
    names(parameters) <- names(params);
    parameters <- paste(names(parameters), parameters, sep=":");
    parameters <- paste(parameters, collapse=";");
    parameters <- paste(celHeader$parameters, parameters, "", sep=";");
    parameters <- gsub(";;", ";", parameters);
    parameters <- gsub(";$", "", parameters);
    celHeader$parameters <- parameters;

    # Create the CEL file
    createCel(pathname, header=celHeader, ..., verbose=less(verbose));
    verbose && exit(verbose);
  } 

  verbose && enter(verbose, "Defining chip-effect file");
  verbose && cat(verbose, "Pathname: ", pathname);
  res <- newInstance(static, pathname);
  verbose && exit(verbose);

  res;
}, static=TRUE)



setMethodS3("readUnits", "ChipEffectFile", function(this, units=NULL, cdf=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- digest(list(units=units, cdf=cdf, ...));
  res <- this$.readUnitsCache[[key]];
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "readUnits.ChipEffectFile(): Returning cached data");
    return(res);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve the data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(cdf)) {
    cdf <- getCellIndices(this, units=units, verbose=less(verbose));
  }

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  res <- NextMethod("readUnits", this, cdf=cdf, ..., verbose=less(verbose));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "readUnits.ChipEffectFile(): Updating cache");
  this$.readUnitsCache <- list();
  this$.readUnitsCache[[key]] <- res;

  res;
})


setMethodS3("getCellIndices", "ChipEffectFile", function(this, ...) {
  getCellIndices(getCdf(this), ...);
})


setMethodS3("updateUnits", "ChipEffectFile", function(this, units=NULL, cdf=NULL, data, ...) {
  if (is.null(cdf))
    cdf <- getCellIndices(this, units=units);

  # Note that the actually call to the encoding is done in updateUnits()
  # of the superclass.
  NextMethod("updateUnits", this, cdf=cdf, data=data, ...);
}, protected=TRUE);


setMethodS3("findUnitsTodo", "ChipEffectFile", function(this, units=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Identifying non-fitted units in chip-effect file");

  verbose && enter(verbose, "Identifying CDF units");
  verbose && cat(verbose, "Pathname: ", getPathname(this));
  verbose && enter(verbose, "Reading CDF cell indices");
  idxs <- getCellIndices(this, units=units, verbose=less(verbose));
  verbose && exit(verbose);
  verbose && enter(verbose, "Extracting first CDF block for each unit");
  idxs <- applyCdfGroups(idxs, .subset2, 1);
  verbose && exit(verbose);
  idxs <- unlist(idxs, use.names=FALSE);
  verbose && exit(verbose);

  # Read one cell from each unit
  verbose && enter(verbose, "Reading data for these cells");
  value <- readCel(getPathname(this), indices=idxs, readIntensities=FALSE, 
                   readStdvs=TRUE, readPixels=FALSE)$stdvs;
  verbose && str(verbose, value);
  verbose && exit(verbose);

  # Identify units for which the stdvs <= 0.
  value <- which(value <= 0);
  if (!is.null(units))
    value <- units[value];
  verbose && str(verbose, value);

  verbose && exit(verbose);

  value;
})


############################################################################
# HISTORY:
# 2006-09-11
# o Great! Using the specially designed CDFs and CEL files to store 
#   estimates is much faster and smaller than using the originally 
#   structured CDF and CEL files.  Now storing the estimates takes a much
#   smaller part of the overall fitting algorithm.
# 2006-09-10
# o Starting to make use of specially design CDFs and CEL files for storing
#   chip effects.  This make getFirstCellIndices() obsolete.
# o Added createParamCdf().
# 2006-08-26
# o Created.  Have to store chip-effect estimates too.  Currently we use
#   the existing CEL/CDF structure for this, but those are unnecessarily
#   large for this.  Later this will be done in special CEL files with a
#   custom CDF file (possible virtual).  This requires that affxparser can
#   create empty CEL files from scratch, which is on the to-do list.
############################################################################
