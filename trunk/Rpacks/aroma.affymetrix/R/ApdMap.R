###########################################################################/**
# @RdocClass ApdMap
#
# @title "The ApdMap class"
#
# \description{
#  @classhierarchy
#
#  An ApdMap object represents an APD read and write map.  
#  Probe indices in an APD data file may be reorder, for instance, such
#  that elements are read in contiguous blocks instead of in an randomized 
#  order which is the typical order in CEL files.  This will speed up reading.
#  An ApdMap object keeps track of the index maps for reading and writing
#  such files.
# }
# 
# @synopsis
#
# \arguments{
#   \item{chipType}{A @character string specifying the chip type.}
#   \item{mapType}{A @character string specifying the map type.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# \seealso{
#  @see "aroma.apd::readApdMap".
# }
#
# @author
#*/###########################################################################
setConstructorS3("ApdMap", function(chipType=NULL, mapType=chipType, ...) {
  extend(Object(), "ApdMap", 
    .chipType = chipType,
    .mapType = mapType,
    .readMap = NULL
  )
})

setMethodS3("as.character", "ApdMap", function(this, ...) {
  s <- paste(class(this)[1], ":", sep="");
  s <- paste(s, " Map type: ", getMapType(this), ".", sep="");
  s;
})


setMethodS3("equals", "ApdMap", function(this, other, ...) {
  if (!identical(this$.chipType, other$.chipType))
    return(FALSE);
  if (!identical(this$.mapType, other$.mapType))
    return(FALSE);
  TRUE;
})


setMethodS3("getChipType", "ApdMap", function(this, ...) {
  this$.chipType;
})


setMethodS3("getMapType", "ApdMap", function(this, ...) {
  this$.mapType;
})


setMethodS3("findApdMap", "ApdMap", function(this, ...) {
  mapType <- this$.mapType;
  mapFile <- findApdMap(mapType);
  mapFile;
})


setMethodS3("getReadMap", "ApdMap", function(this, ...) {
  readMap <- this$.readMap;
  if (is.null(readMap)) {
    mapType <- this$.mapType;
    mapFile <- findApdMap(mapType);
    if (is.null(mapFile))
      throw("Could not find APD map file for map type: ", mapType);
    readMap <- readApdMap(mapFile)$map;
    this$.readMap <- readMap;
  }
  readMap;  
})

setMethodS3("getWriteMap", "ApdMap", function(this, ...) {
  invertMap(getReadMap(this));
})


setMethodS3("read", "ApdMap", function(static, filename, path=NULL, ...) {
  # Argument 'filename' and 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);

  map <- readApdHeader(pathname);
  chipType <- map$chipType;
  mapType <- basename(pathname);  
  mapType <- gsub("[.]apm$", "", mapType);

  map <- newInstance(static, chipType=chipType, mapType=mapType, ...);
  readMap <- getReadMap(map);
  map;
}, static=TRUE)



setMethodS3("write", "ApdMap", function(this, path=NULL, overwrite=FALSE, ...) {
  mapType <- this$.mapType;
  mapFile <- paste(mapType, ".apm", sep="");
  pathname <- Arguments$getWritablePathname(mapFile, path=path, 
                                                   mustNotExist=!overwrite);
  readMap <- this$.readMap;
  if (is.null(readMap)) {
    mapFile <- findApdMap(mapType);
    if (is.null(mapFile))
      throw("Could not find APD map file for map type: ", mapType);
  }

  # Write read map to file
  writeApdMap(pathname, map=readMap, chipType=getChipType(this));

  invisible(pathname);
})


setMethodS3("fromMapType", "ApdMap", function(static, mapType, ...) {
  # Argument 'mapType':
  mapFile <- findApdMap(mapType);
  if (is.null(mapFile))
    throw("Could not find APD map file for map type: ", mapType);
  chipType <- readApdHeader(mapFile)$chipType;
  newInstance(static, mapType=mapType, chipType=chipType, ...);
}, static=TRUE)


setMethodS3("fromCdf", "ApdMap", function(static, chipType, label=NULL, force=FALSE, ...) {
  # Argument 'label':
  label <- Arguments$getCharacter(label);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'chipType':
  cdfFile <- findCdf(chipType);
  if (is.null(cdfFile)) {
    throw("Could not find CDF file for chip type: ", chipType);
  }

  # Generate the pathname for the map file
  mapType <- chipType;
  if (!is.null(label))
    mapType <- paste(mapType, "-", label, sep="");

  apdMap <- ApdMap(chipType=chipType, mapType=mapType);

  # Map file already exists?
  mapFile <- findApdMap(mapType);
  if (!is.null(mapFile) && !force) {
    return(apdMap);
  }

  # Generate optimal read map
  writeMap <- readCdfUnitsWriteMap(cdfFile, ...);   # in [1,n]
  apdMap$.readMap <- invertMap(writeMap);
  
  apdMap;
}, static=TRUE);




############################################################################
# HISTORY:
# 2006-05-29
# o Added equals().
# 2006-05-15
# o Created.
############################################################################
