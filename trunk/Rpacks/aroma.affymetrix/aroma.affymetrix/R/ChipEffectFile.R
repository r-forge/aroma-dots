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
#   \item{probeModel}{The specific type of model, e.g. \code{"pm"}.}
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
setConstructorS3("ChipEffectFile", function(..., probeModel=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probeModel':
  probeModel <- match.arg(probeModel);

  this <- extend(ParameterCelFile(...), "ChipEffectFile",
    "cached:.firstCells" = NULL,
    probeModel = probeModel
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


setMethodS3("clearCache", "ChipEffectFile", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".firstCells")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)


setMethodS3("as.character", "ChipEffectFile", function(this, ...) {
  s <- NextMethod(generic="as.character", object=this, ...);
  params <- paste(getParametersAsString(this), collapse=", ");
  s <- c(s, sprintf("Parameters: (%s)", params));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getParameters", "ChipEffectFile", function(this, ...) {
  params <- list(
    probeModel = this$probeModel
  );
  params;
})

setMethodS3("getParametersAsString", "ChipEffectFile", function(this, ...) {
  params <- getParameters(this);
  params <- trim(capture.output(str(params)))[-1];
  params <- gsub("^[$][ ]*", "", params);
  params <- gsub(" [ ]*", " ", params);
  params <- gsub("[ ]*:", ":", params);
  params;
}, private=TRUE)


setMethodS3("createParamCdf", "ChipEffectFile", function(static, sourceCdf, ..., verbose=FALSE) {
  # Argument 'verbose': 
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Creating CDF for chip effects");
  verbose && cat(verbose, "Source chip type: ", getChipType(sourceCdf));
  verbose && cat(verbose, "Source CDF: ", getPathname(sourceCdf));

  # Search for existing monocell CDF
  for (sep in c(",", "-")) {
    chipType <- paste(getChipType(sourceCdf), "monocell", sep=sep);
    verbose && cat(verbose, "Looking for chip type: ", chipType);
    pathname <- AffymetrixCdfFile$findByChipType(chipType);
    if (!is.null(pathname)) {
      verbose && cat(verbose, "Found: ", pathname);
      break;
    }
  }
  # Warn about deprecated filname <chipType>-monocell.
  if (!is.null(pathname) && (sep == "-")) {
    msg <- paste("Deprecated filename of monocell CDF detected (uses dash instead of comma): ", pathname);
    warning(msg);
    verbose && cat(verbose, msg);
    verbose && enter(verbose, "Renaming (old-style) monocell CDF");
    verbose && cat(verbose, "Source: ", pathname);
    dest <- gsub("-monocell[.]", ",monocell.", pathname);
    verbose && cat(verbose, "Destination: ", dest);
    res <- file.rename(pathname, dest);
    if (!res)
      throw("Failed to rename monocell CDF file: ", pathname, " -> ", dest);
    pathname <- dest;
    verbose && exit(verbose, msg);
  }

  if (is.null(pathname)) {
    verbose && cat(verbose, "Pathname: Not found!");
    verbose && cat(verbose, "Will create CDF for the chip-effect files from the original CDF. NOTE: This will take several minutes or more!");
    verbose && enter(verbose, "Creating CDF");
    cdf <- createMonoCell(sourceCdf, verbose=less(verbose));
    verbose && exit(verbose);
  } else {
    verbose && cat(verbose, "Pathname: ", pathname);
    cdf <- AffymetrixCdfFile$fromFile(pathname);
  }
  verbose && exit(verbose);

  cdf;
}, static=TRUE, private=TRUE)


setMethodS3("fromDataFile", "ChipEffectFile", function(static, df=NULL, filename=sprintf("%s,chipEffects.cel", getFullName(df)), path, name=getName(df), cdf=NULL, ..., verbose=FALSE) {
  # Argument 'df':
  if (!is.null(df)) {
    if (!inherits(df, "AffymetrixCelFile"))
      throw("Argument 'df' is not an AffymetrixCelFile: ", class(df)[1]);
  }

  # Argument 'cdf':
  if (is.null(cdf)) {
    if (is.null(df))
      throw("Either argument 'df' or 'cdf' must specified.");
  } else {
    if (!inherits(cdf, "AffymetrixCdfFile"))
      throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }


  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  pathname <- Arguments$getWritablePathname(filename, path=path);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Backward compatibility patch for now. Before chip effect files
  # only carried on the sample name, but not the tags. If such a 
  # file is detected, it is renamed. 
  # This should be removed in future versions. /HB 2007-01-10
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!isFile(pathname)) {
    filenameOld <- sprintf("%s,chipEffects.cel", getName(df));
    pathnameOld <- Arguments$getWritablePathname(filenameOld, path=path);
    if (isFile(pathnameOld)) {
      verbose && enter(verbose, "Renaming chip-effect file using old style file format");
      verbose && cat(verbose, "Source: ", pathnameOld);
      verbose && cat(verbose, "Destination: ", pathname);
      file.rename(pathnameOld, pathname);
      if (!isFile(pathname)) {
        throw("Failed to rename chip-effect from old to new filename format: ", pathnameOld);
      }
      verbose && exit(verbose);
    }
    rm(filenameOld, pathnameOld);
  }


  if (!isFile(pathname)) {
    verbose && enter(verbose, "Creating chip-effect file");
    verbose && cat(verbose, "Pathname: ", pathname);

    # Get CDF for chip effects
    if (is.null(cdf))
      cdf <- createParamCdf(static, getCdf(df), verbose=less(verbose));
  
    # Get CDF header
    cdfHeader <- getHeader(cdf);

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

##    # Fill with negative values
##    nbrOfProbes <- celHeader$total;
##    updateCel(pathname, indices=1:nbrOfProbes, intensities=rep(-1,nbrOfProbes), verbose=less(verbose));

    verbose && exit(verbose);
  } 

  verbose && enter(verbose, "Defining chip-effect file");
  verbose && cat(verbose, "Pathname: ", pathname);
  res <- newInstance(static, pathname);
  verbose && exit(verbose);

  res;
}, static=TRUE, private=TRUE)



setMethodS3("readUnits", "ChipEffectFile", function(this, units=NULL, cdf=NULL, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # Check for cached data
  key <- list(method="readUnits", class=class(this)[1], 
              pathname=getPathname(this),
              cdf=cdf, units=units, ...);
  id <- digest(key);
  res <- this$.readUnitsCache[[id]];
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
  res <- NextMethod("readUnits", this, cdf=cdf, ..., force=force, verbose=less(verbose));

  # Store read units in cache?
  if (cache) {
    verbose && cat(verbose, "readUnits.ChipEffectFile(): Updating cache");
    this$.readUnitsCache <- list();
    this$.readUnitsCache[[id]] <- res;
  }

  res;
})


setMethodS3("getCellIndices", "ChipEffectFile", function(this, ..., .cache=TRUE) {
  getCellIndices(getCdf(this), ...);
})


setMethodS3("updateUnits", "ChipEffectFile", function(this, units=NULL, cdf=NULL, data, ...) {
  if (is.null(cdf))
    cdf <- getCellIndices(this, units=units);

  # Note that the actually call to the encoding is done in updateUnits()
  # of the superclass.
  NextMethod("updateUnits", this, cdf=cdf, data=data, ...);
}, private=TRUE);



setMethodS3("findUnitsTodo", "ChipEffectFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Identifying non-fitted units in chip-effect file");

  verbose && cat(verbose, "Pathname: ", getPathname(this));


  idxs <- NULL;
  if (is.null(units)) {
    # Look up chip-type and parameter specific but data set independent data
    cdf <- getCdf(this);
    chipType <- getChipType(cdf);
    key <- list(method="findUnitsTodo", class=class(this)[1], 
                chipType=chipType, params=getParameters(this));
    dirs <- c("aroma.affymetrix", chipType);
    if (!force) {
      idxs <- loadCache(key, dirs=dirs);
      if (!is.null(idxs))
        verbose && cat(verbose, "Found indices cached on file");
    }
  }

  if (is.null(idxs)) {
    verbose && enter(verbose, "Identifying CDF units");
  
    verbose && enter(verbose, "Reading CDF cell indices");
    idxs <- getCellIndices(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  
    verbose && enter(verbose, "Extracting first CDF block for each unit");
    idxs <- applyCdfGroups(idxs, .subset2, 1);
    verbose && exit(verbose);
  
    idxs <- unlist(idxs, use.names=FALSE);

    if (is.null(units)) {
      verbose && enter(verbose, "Saving to file cache");
      saveCache(idxs, key=key, dirs=dirs);
      verbose && exit(verbose);
    }

    verbose && exit(verbose);
  }


  # Read one cell from each unit
  verbose && enter(verbose, "Reading data for these ", length(idxs), " cells");
  value <- readCel(getPathname(this), indices=idxs, readIntensities=FALSE, 
                   readStdvs=TRUE, readPixels=FALSE)$stdvs;
  verbose && exit(verbose);


  # Identify units for which the stdvs <= 0.
  value <- which(value <= 0);
  if (!is.null(units))
    value <- units[value];
  verbose && cat(verbose, "Looking for stdvs <= 0 indicating non-estimated units:");
  verbose && str(verbose, value);

  verbose && exit(verbose);

  value;
})



setMethodS3("getCellMap", "ChipEffectFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (inherits(units, "ChipEffectFileCellMap")) {
    return(units);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving unit-to-cell map");

  # Is 'units' already a CDF list?
  if (is.list(units)) {
    # No fancy validation for now.
    cells <- units;
    cdf <- getCdf(this);
    units <- indexOf(cdf, names=names(units));
    if (any(is.na(units))) {
      throw("Argument 'units' is of unknown structure.");
    }
    verbose && enter(verbose, "Argument 'cells' is already a CDF cell-index structure");
  } else {
    verbose && enter(verbose, "Retrieving cell indices for specified units");
    # Get the cells to read
    cells <- getCellIndices(this, units=units, force=force, verbose=less(verbose));
  }

  unitNames <- names(cells);
  unitSizes <- unlist(lapply(cells, length), use.names=FALSE);
  cells <- unlist(cells, use.names=FALSE);
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Creating return data frame");
  uUnitSizes <- unique(unitSizes);
  if (is.null(units)) {
    cdf <- getCdf(this);
    units <- seq(length=nbrOfUnits(cdf));
  }
  units <- rep(units, each=unitSizes);

  # The following is too slow:
  #  groups <- sapply(unitSizes, FUN=function(n) seq(length=n));

  # Instead, updated size by size
  groups <- matrix(NA, nrow=max(uUnitSizes), ncol=length(unitNames));
  for (size in uUnitSizes) {
    cc <- which(unitSizes == size);
    seq <- seq(length=size);
    groups[seq,cc] <- seq;
  }
  groups <- groups[!is.na(groups)];
  map <- data.frame(unit=units, group=groups, cell=cells);
  verbose && exit(verbose);

  verbose && exit(verbose);

  class(map) <- c("ChipEffectFileCellMap", class(map));

  map;
}, private=TRUE)



setMethodS3("getDataFlat", "ChipEffectFile", function(this, units=NULL, fields=c("theta", "sdTheta", "outliers"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving data as a flat data frame");

  # Get unit-to-cell map
  suppressWarnings({
    map <- getCellMap(this, units=units, ..., verbose=less(verbose));
  })

  verbose && enter(verbose, "Reading data fields");
  celFields <- c(theta="intensities", sdTheta="stdvs", outliers="pixels");
  suppressWarnings({
    data <- getData(this, indices=map[,"cell"], fields=celFields[fields]);
  })
  rownames(data) <- seq(length=nrow(data));  # Work around?!? /HB 2006-11-28

  # Decode
  names <- colnames(data);
  names <- gsub("intensities", "theta", names);
  names <- gsub("stdvs", "sdTheta", names);
  names <- gsub("pixels", "outliers", names);
  colnames(data) <- names;
  verbose && str(verbose, data);
  if ("outliers" %in% names) {
    data[,"outliers"] <- as.logical(-data[,"outliers"]);
  }
  verbose && exit(verbose);

  len <- sapply(data, FUN=length);
  keep <- (len == nrow(map));
  data <- data[keep];
  data <- as.data.frame(data);

  data <- cbind(map, data);

  verbose && exit(verbose);

  data;
}, private=TRUE)



setMethodS3("updateDataFlat", "ChipEffectFile", function(this, data, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'data':
  names <- colnames(data);  
  namesStr <- paste(names, collapse=", ");
  if (!"cell" %in% names)
    throw("Argument 'data' must contain a column 'cell'": namesStr);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose2 <- -as.integer(verbose)-2;

  verbose && enter(verbose, "Storing flat data to file");

  # Encode
  names <- gsub("theta", "intensities", names);
  names <- gsub("sdTheta", "stdvs", names);
  names <- gsub("outliers", "pixels", names);
  colnames(data) <- names;
  if ("pixels" %in% names) {
    data[,"pixels"] <- -as.integer(data[,"pixels"]);
  }

  verbose && enter(verbose, "Updating file");
  indices <- data[,"cell"];
  keep <- (names %in% c("intensities", "stdvs", "pixels"));
  data <- data[,keep];
  pathname <- getPathname(this);
  updateCel(pathname, indices=indices, data, verbose=verbose2);
  verbose && exit(verbose);

  verbose && exit(verbose);
  invisible(data);
}, private=TRUE)



############################################################################
# HISTORY:
# 2007-02-09
# o Updated the file cache sub directory.
# 2007-01-10
# o Now fromDataFile() looks for chip effect files named using the "sample
#   name" only (not tags) file format, and renames it to the full name
#   format.  This is a "patch" so we don't have to reestimate old data sets.
# 2007-01-09
# o Now fromDataFile() generates a file with the full name (name + tags) of
#   the input file and not just the name.
# 2007-01-07
# o TO DO: Speed up getCellMap(), e.g. using more clever file caching.
# 2007-01-06
# o findUnitsTodo() cache cell-index vector to file if all units are to
#   be scanned.
# 2007-01-05
# o Removed getSampleNames().
# 2007-10-02
# o TO DO: Static fromDataFile() does not really need argument 'df'.
# 2006-12-02
# o BUG FIX: getCellMap(..., units=NULL) did not work.
# 2006-11-28
# o Added trial version of updateDataFlat(). Seems to work. Will speed
#   up a few things.
# o Added trial versions of getCellMap() and getDataFlat().
# 2006-11-28
# o Added argument 'cache' to readUnits() to specify if the result should
#   be cached or not.
# o BUG FIX: The caching mechanism of readUnits() of ChipEffectFile was
#   not sensitive to the class of the ChipEffectFile object.  Added a
#   'class' element to the cache key.  Then each subclass must override
#   this method too.
# 2006-10-06
# o Now chip effect files use filename tag 'chipEffects' with a comma in
#   front (instead of suffix "-chipEffects').  This way getName() will
#   return the sample name without any extra endings.
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
