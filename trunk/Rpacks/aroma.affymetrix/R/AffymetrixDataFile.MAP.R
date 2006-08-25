# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 
#  O P T I M I Z E D   R E A D   &   W R I T E   M A P 
# 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
###########################################################################/**
# @set "class=AffymetrixDataFile"
# @RdocMethod getMapType
#
# @title "Gets the map type of the data file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getMapType", "AffymetrixDataFile", abstract=TRUE, protected=TRUE);



###########################################################################/**
# @RdocMethod getReadMap
#
# @title "Gets the probe read map for a data file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{inverse}{If @FALSE, the read map is returned, otherwise the write map.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "setApdMap".
#   @seemethod "hasApdMap".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getReadMap", "AffymetrixDataFile", function(this, inverse=FALSE, ...) {
  apdMap <- this$.apdMap;
  if (is.null(apdMap))
    return(NULL);

  readMap <- getReadMap(apdMap);

  # Inverse map?
  if (inverse)
    readMap <- invertMap(readMap);

  readMap;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod setApdMap
#
# @title "Sets the probe read map for a data file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{apdMap}{An @see "ApdMap" object or @NULL.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "getReadMap".
#   @seemethod "hasApdMap".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setApdMap", "AffymetrixDataFile", function(this, apdMap=NULL, ...) {
  # Argument 'apdMap':
  if (is.null(apdMap)) {
  } else if (!inherits(apdMap, "ApdMap")) {
    throw("Argument 'apdMap' is of unknown type: ", class(apdMap)[1]);
  }

  this$.apdMap <- apdMap;
  invisible(this);
}, protected=TRUE)


setMethodS3("getApdMap", "AffymetrixDataFile", function(this, ...) {
  apdMap <- this$.apdMap;
  if (is.null(apdMap)) {
    mapType <- getMapType(this);
    if (!is.null(mapType)) {
      apdMap <- ApdMap$fromMapType(mapType, ...);
      this$.apdMap <- apdMap;
    }
  }
  apdMap;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod hasApdMap
#
# @title "Checks if a data file has a probe map"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE if a read map exists, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seemethod "getReadMap".
#   @seemethod "setApdMap".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("hasApdMap", "AffymetrixDataFile", function(this, ...) {
  !is.null(this$.apdMap);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod readCdfUnitsWriteMap
#
# @title "Read the probe map for a subset of units"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{units}{A @numeric @vector of indices of units to be read. 
#    If @NULL, all units are considered.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer @vector which is a \emph{write} map.
# }
#
# @author
#
# \seealso{
#   Internally @see "affxparser::readCdfUnitsWriteMap" is used.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("readCdfUnitsWriteMap", "AffymetrixDataFile", function(this, units=NULL, ...) {
  cdfFile <- findCdf(this);
  writeMap <- readCdfUnitsWriteMap(cdfFile, units=units);   # in [1,n]

  writeMap;
}, protected=TRUE)



setMethodS3("createOptimalReadMap", "AffymetrixDataFile", function(this, path=NULL, label=NULL, force=FALSE, ...) {
  # Argument 'path':
  if (!is.null(path)) {
    path <- Arguments$getReadablePath(path, mustExist=TRUE);
  }

  # Argument 'label':
  label <- Arguments$getCharacter(label);

  # Argument 'force':
  force <- Arguments$getLogical(force);


  # Generate the pathname for the map file
  chipType <- getChipType(this);
  mapType <- chipType;
  if (!is.null(label))
    mapType <- paste(mapType, "-", label, sep="");
  mapFile <- paste(mapType, ".apm", sep="");
  pathname <- filePath(path, mapFile);

  # Map file already exists?
  if (file.exists(pathname) && !force) {
    readMap <- readApdMap(pathname);
    return(mapType);
  }

  # Generate optimal read map
  cdfFile <- findCdf(this);
  writeMap <- readCdfUnitsWriteMap(cdfFile, ...);   # in [1,n]
  readMap <- invertMap(writeMap);

  # Write read map to file
  writeApdMap(pathname, map=readMap);

  mapType;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2006-05-15
# o Renamed getMap() to getReadMap() and same for the other map functions.
#   Removed mapOn(), mapOff() etc.
# o Extracted from AffymetrixDataFile.R.
# 2006-04-09
# o Added abstract method getMapType().
# 2006-03-04
# o Remapping of indices works only with integers, which means that we can
#   only index 256^4 = 4.2*10^9 cells.
# o Added get- and setMap() and hasMap().
# o Now getProbeSignals() and readIntensities() accepts argument 'map'.
############################################################################
