setMethodS3("getUnitGroupCellMap", "AffymetrixCdfFile", function(this, units=NULL, retNames=FALSE, force=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  flattenCellIndices <- function(cells, ..., retNames=FALSE, verbose=FALSE) { 
    # Returning indices or names?
    if (!retNames) {
      verbose && enter(verbose, "Renaming group names to group indices");
      cells <- base::lapply(cells, FUN=function(unit) {
        groups <- .subset2(unit, 1);
        names(groups) <- seq_len(length(groups));
        list(groups=groups);
      });

      verbose && print(verbose, gc);
      verbose && exit(verbose);
    }

    # Flatten cell data 
    verbose && enter(verbose, "Flattening cell data");
    cells <- unlist(cells, use.names=TRUE);
    names <- names(cells);
    names(cells) <- NULL;
    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);

    verbose && enter(verbose, "Extract unit and group names");
    # Do some tricks to clean up the names 
    names <- gsub("([.]groups|indices*)", "", names);
    pattern <- "^(.*)[.](.*)[.](.*)$";

    units <- gsub(pattern, "\\1", names);
    groups <- gsub(pattern, "\\2", names);
    rm(pattern, names); # Not needed anymore
    verbose && exit(verbose);

    if (!retNames) {
      verbose && enter(verbose, "Converting to indices");
      units <- as.integer(units);
      groups <- as.integer(groups);
    }
    verbose && exit(verbose);
  
    # Return data 
    if (retNames) {
      map <- data.frame(unit=units, group=groups, cell=cells);
    } else {
      map <- matrix(c(units, groups, cells), ncol=3, 
                         dimnames=list(NULL, c("unit", "group", "cell")));
    }

    map;
  } # flattenCellIndices()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'retNames':
  retNames <- Arguments$getLogical(retNames);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Getting (unit, group, cell) map");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this); 

  # Look for results in file cache 
  verbose && enter(verbose, "Checking cache");
  key <- list(method="getUnitGroupCellMap", class=class(this)[1], 
                   chipType=chipType, units=units, retNames=retNames, ...); 
  dirs <- c("aroma.affymetrix", chipType); 
  map <- NULL;
  if (!force) 
     map <- loadCache(key=key, dirs=dirs); 
  if (is.null(map)) {
    verbose && exit(verbose, suffix="...miss");
  } else {
    verbose && printf(verbose, "RAM: %.2fMB\n", object.size(map)/1024^2);
    verbose && exit(verbose, suffix="...hit");
  }

  # Not in cache?
  if (is.null(map)) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Retrieve matrix of (unit, group, cell) indices
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    cells <- getCellIndices(this, units=units, ..., verbose=less(verbose));
    nbrOfUnits <- length(cells);
    verbose && printf(verbose, "Read %d units\n", nbrOfUnits);
  
    if (!retNames) {
      # Convert unit names to unit indices
      if (is.null(units))
        units <- seq_len(nbrOfUnits);
      names(cells) <- units;
  
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
    }
  
    verbose && enter(verbose, "Flattening cell indices to create cell map");
    map <- flattenCellIndices(cells, retNames=retNames, verbose=less(verbose)); 
    if (!retNames) {
      map <- as.matrix(map);
      rownames(map) <- NULL;
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
    }
    verbose && str(verbose, map);
    verbose && exit(verbose);
  
  
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Save to cache
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Save only results > 500kB
    if (object.size(map) > 500e3) {
      saveCache(map, key=key, dirs=dirs); 
      verbose && cat(verbose, "Saved to file cache");
    }
  } # if (is.null(map))

  # Extract subset of units
##  if (!is.null(units)) {
##    mapUnits <- map[,"unit"];
##    if (retNames) {
##      # Convert unit names in the map to unit indices
##      allUnitNames <- getUnitNames(this);
##      mapUnits <- match(mapUnits, allUnitNames);
##    }
##    keep <- (mapUnits %in% units);
##    map <- map[keep,,drop=FALSE];
##  }
  
  verbose && printf(verbose, "RAM: %.2fMB\n", object.size(map)/1024^2);
  verbose && exit(verbose);

  map; 
}, protected=TRUE)  # getUnitGroupCellMap()



############################################################################
# HISTORY:
# 2007-03-08
# o Note, getUnitGroupCellMap() is just a test function.
# o Created.
############################################################################
