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
# }
#
#*/###########################################################################
setConstructorS3("ChipEffectFile", function(..., model=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  extend(ParameterCelFile(...), "ChipEffectFile",
    "cached:.firstCells" = NULL,
    model = model
  )
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
    cdf <- createMonoCell(sourceCdf, verbose=verbose);
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

  pathname <- Arguments$getWritablePathname(filename, path=path);
  if (!isFile(pathname)) {
    # Get CDF for chip effects
    if (is.null(cdf))
      cdf <- createParamCdf(static, getCdf(df), verbose=verbose);
  
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

#     celHeader0 <- celHeader;
#     celHeader$rows <- cdfHeader$rows;
#     celHeader$cols <- cdfHeader$cols;
#     celHeader$total <- cdfHeader$total;
#     celHeader$chiptype <- cdfHeader$chiptype;

    # Create the CEL file
    createCel(pathname, header=celHeader, ..., verbose=verbose);
  }

  newInstance(static, pathname);
}, static=TRUE)



setMethodS3("encodeUnitGroup", "ChipEffectFile", function(static, groupData, ...) {
  theta <- .subset2(groupData, "theta");
  stdvs <- .subset2(groupData, "sdTheta");
  outliers <- .subset2(groupData, "thetaOutliers");
  pixels <- -as.integer(outliers);

  list(intensities=theta, stdvs=stdvs, pixels=pixels);
}, static=TRUE, protected=TRUE)


setMethodS3("decodeUnitGroup", "ChipEffectFile", function(static, groupData, ...) {
  res <- list();
  if (!is.null(groupData$intensities))
    res$theta <- groupData$intensities;
  if (!is.null(groupData$stdvs))
    res$sdTheta <- groupData$stdvs;
  if (!is.null(groupData$pixels))
    res$thetaOutliers <- as.logical(-groupData$pixels);
  res;
}, static=TRUE, protected=TRUE)


setMethodS3("readUnits", "ChipEffectFile", function(this, units=NULL, cdf=NULL, ...) {
  if (is.null(cdf)) {
    # Use only the first cell in each unit group.
    cdf <- readCdfCellIndices(getPathname(getCdf(this)), units=units);
  }

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  NextMethod("readUnits", this, cdf=cdf, ...);
})


setMethodS3("updateUnits", "ChipEffectFile", function(this, units=NULL, cdf=NULL, data, ...) {
  if (is.null(cdf)) {
    # Use only the first cell in each unit group.
    cdf <- readCdfCellIndices(getPathname(getCdf(this)), units=units);
  }

  # Note that the actually call to the encoding is done in updateUnits()
  # of the superclass.
  NextMethod("updateUnits", this, cdf=cdf, data=data, ...);
}, protected=TRUE);


setMethodS3("findUnitsTodo", "ChipEffectFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  idxs <- getCellIndices(getCdf(this), units=NULL);
  idxs <- applyCdfGroups(idxs, .subset2, 1);
  idxs <- unlist(idxs, use.names=FALSE);

  # Read one cell from each unit
  verbose && enter(verbose, "Reading data for these cells");
  value <- readCel(getPathname(this), indices=idxs, readIntensities=FALSE, 
                   readStdvs=TRUE, readPixels=FALSE)$stdvs;
  verbose && str(verbose, value);
  verbose && exit(verbose);

  # Identify units for which the stdvs <= 0.
  which(value <= 0);
})


############################################################################
# HISTORY:
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
