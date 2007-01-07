###########################################################################/**
# @RdocClass PlatformDesign
#
# @title "The PlatformDesign class"
#
# \description{
#  @classhierarchy
#
#  The PlatformDesign class provides a memory-efficient interface to a
#  so called platform-design package.  Its method can access a subset of the
#  data of such packages without having to load the package, which saves
#  memory.
# }
# 
# @synopsis
#
# \arguments{
#  \item{cdf}{An @see "AffymetrixCdfFile".}
#  \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("PlatformDesign", function(cdf=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(cdf)) {
    if (!inherits(cdf, "AffymetrixCdfFile")) {
      throw("Argument 'cdf' is not an AffymetrixCdfFile object: ", 
                                                              class(cdf)[1]);
    }
  }

  extend(Object(), "PlatformDesign",
    .cdf = cdf
  )
}, private=TRUE)


setMethodS3("as.character", "PlatformDesign", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  s <- c(s, sprintf("Package name: %s", getPackageName(this)));
  s <- c(s, sprintf("Is loaded: %s", isLoaded(this)));
  if (isInstalled(this)) {
    s <- c(s, sprintf("Version: %s", getDescription(this)$Version));
    s <- c(s, sprintf("Date: %s", getDescription(this)$Created));
    s <- c(s, sprintf("Path: %s", getPath(this)));
  } else {
    s <- c(s, "WARNING: Package is not installed!");
  }
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("isInstalled", "PlatformDesign", function(this, ...) {
  path <- getPath(this);
  (path != "");
})

setMethodS3("isLoaded", "PlatformDesign", function(this, ...) {
  pkgName <- getPackageName(this);
  pattern <- sprintf("^package:%s$", pkgName);
  any(regexpr(pattern, search()) != -1);
})

setMethodS3("load", "PlatformDesign", function(this, ...) {
  pkgName <- getPackageName(this);
  if (!require(pkgName, character.only=TRUE)) {
    throw("Package not loaded: ", pkgName);
  }
})

setMethodS3("unload", "PlatformDesign", function(this, ...) {
  pkgName <- getPackageName(this);
  pattern <- sprintf("^package:%s$", pkgName);
  while (TRUE) {
    pos <- grep(pattern, search());
    if (length(pos) == 0)
      break;
    detach(pos=pos[1]);
  }
})

setMethodS3("getCdf", "PlatformDesign", function(this, ...) {
  this$.cdf;
})

setMethodS3("getChipType", "PlatformDesign", function(this, clean=FALSE, ...) {
  cdf <- getCdf(this);
  chipType <- getChipType(cdf);
  if (clean) {
    chipType <- cleancdfname(chipType, addcdf=FALSE);
  }
  chipType;
})

setMethodS3("getName", "PlatformDesign", function(this, ...) {
  getChipType(this, clean=TRUE);
})

setMethodS3("getPackageName", "PlatformDesign", function(where, ...) {
  # To please R CMD check
  this <- where;
  sprintf("pd%s", getChipType(this, clean=TRUE));
})

setMethodS3("getPath", "PlatformDesign", function(this, ...) {
  system.file(package=getPackageName(this))
})

setMethodS3("getDescription", "PlatformDesign", function(this, ...) {
  packageDescription(getPackageName(this));
})

setMethodS3("getDataFiles", "PlatformDesign", function(this, ...) {
  path <- file.path(getPath(this), "data");
  files <- list.files(path, all.files=TRUE);
  files <- files[(regexpr("([.]|[.][.])$", files) == -1)];
  files <- file.path(path, files);
  files;
}, private=TRUE)

setMethodS3("getDataFile", "PlatformDesign", function(this, pattern, ...) {
  path <- file.path(getPath(this), "data");
  list.files(path=path, pattern=pattern, full.names=TRUE);
}, private=TRUE)

setMethodS3("loadDataFile", "PlatformDesign", function(this, filename, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'filename':
  pathname <- file.path(getPath(this), "data", filename);
  if (!isFile(pathname)) {
    throw("File not found: ", pathname);
  }

  verbose && printf(verbose, "Pathname: %s\n", pathname);
  filesize <- file.info(pathname)$size;
  verbose && printf(verbose, "File size: %.1fMB\n", filesize/1024^2);

  env <- new.env();  
  vars <- base::load(file=pathname, envir=env);
  attr(env, "pathname") <- pathname;
  attr(env, "vars") <- vars;

  env;
}, private=TRUE)


setMethodS3("getAnnotations", "PlatformDesign", function(this, ...) {
  filename <- sprintf("%s.rda", getName(this));
  env <- loadDataFile(this, filename=filename, ...);
  var <- attr(env, "vars");
  env <- env[[var]];
  env;
}, private=TRUE)

setMethodS3("getCrlmmInfo", "PlatformDesign", function(this, ...) {
  filename <- sprintf("%sCrlmmInfo.rda", getName(this));
  env <- loadDataFile(this, filename=filename, ...);
  var <- attr(env, "vars");
  env <- env[[var]];
  env;
}, private=TRUE)


setMethodS3("getPlatformDesignObject", "PlatformDesign", function(this, ...) {
  filename <- sprintf("pd%s.rda", getName(this));
  env <- loadDataFile(this, filename=filename, ...);
  var <- attr(env, "vars");
  env <- env[[var]];
  env;
}, private=TRUE)


setMethodS3("getReferenceQuantiles", "PlatformDesign", function(this, ...) {
  filename <- sprintf("pd%sRef.rda", getName(this));
  env <- loadDataFile(this, filename=filename, ...);
  var <- attr(env, "vars");
  env <- env[[var]];
  env;
})


setMethodS3("getProbeEffects", "PlatformDesign", function(this, ...) {
  filename <- sprintf("pd%sRmaPLM.rda", getName(this));
  env <- loadDataFile(this, filename=filename, ...);
  var <- attr(env, "vars");
  env <- env[[var]];
  env;
})

setMethodS3("getCachePath", "PlatformDesign", function(this, fields, ...) {
  # Get the path to the file cache for this platform design
  path <- filePath(getCachePath(), getName(this));
  path;
}, private=TRUE)


setMethodS3("getCachedFields", "PlatformDesign", function(this, fields, ...) {
  path <- getCachePath(this);
  files <- list.files(path=path);
  files <- gsub("[.]rda$", "", files);
  files;
}, private=TRUE)


setMethodS3("clearFileCache", "PlatformDesign", function(this, fields, ...) {
  path <- getCachePath(this);
  files <- list.files(path=path);
  for (file in files) {
    pathname <- file.path(path, file);
    file.remove(pathname);
  }
  invisible(files);
})


setMethodS3("getFeatureInfo", "PlatformDesign", function(this, fields, subset=NULL, drop=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Get the path to the file cache for this platform design
  path <- getCachePath(this);
  mkdirs(path);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check what fields are in the file cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filenames <- sprintf("%s.rda", fields);
  pathnames <- file.path(path, filenames);
  exists <- file.exists(pathnames);

  if (all(exists)) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Load all data from file cache
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    env <- new.env();
    for (kk in seq(along=fields)) {
      field <- fields[kk];
      value <- loadObject(pathnames[kk]);
      if (!is.null(subset))
        value <- value[subset];
      assign(field, value, envir=env);
      rm(value);
    }
  } else {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Load all data from platform-design package and store all in file cache
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Load missing data
    env <- getPlatformDesignObject(this, verbose=verbose);
    env <- featureInfo(env);

    verbose && enter(verbose, "Store all new fields in file cache");
    # Store new fields to the file cache 
    newFields <- setdiff(ls(envir=env), fields[exists]);
    for (kk in seq(along=newFields)) {
      field <- newFields[kk];
      filename <- sprintf("%s.rda", field);
      pathname <- file.path(path, filename);
      verbose && enter(verbose, "Saving field '", field, "'");
      verbose && cat(verbose, "Pathname: ", pathname);
      value <- env[[field]];
      saveObject(value, file=pathname);
      rm(value);
      verbose && exit(verbose);
    }
    verbose && exit(verbose);

    # Remove unwanted fields
    rm(list=setdiff(ls(envir=env), fields), envir=env);
    gc(); # Here it is ok to call gc() because it should happen too often!

    # Extract subset
    if (!is.null(subset)) {
      for (kk in seq(along=fields)) {
        field <- fields[kk];
        env[[field]] <- env[[field]][subset];
      }
    }
    gc();
  }

  if (drop && length(fields) == 1) {
    env <- env[[fields]];
  }

  env;
})

setMethodS3("getUnitsOnChromosome", "PlatformDesign", function(this, chromosome, ...) {
  str <- paste("chr", chromosome, sep="");
  which(getAnnotations(this)$Chromosome == str);
})

setMethodS3("getCellIndices", "PlatformDesign", function(this, ...) {
  # The feature index as indexed by the CDF. Totally redundant with (X,Y).
  as.integer(getFeatureInfo(this, "order_index", ...));
})

setMethodS3("getReadMap", "PlatformDesign", function(this, ...) {
  map <- this$.readMap;
  if (is.null(map)) {
    map <- affxparser::invertMap(getCellIndices(this));
    this$.readMap <- map;
  }
  map;
}, private=TRUE)

setMethodS3("getWriteMap", "PlatformDesign", function(this, ...) {
  affxparser::invertMap(getReadMap(this));
}, private=TRUE)

setMethodS3("getX", "PlatformDesign", function(this, ...) {
  as.integer(getFeatureInfo(this, "X", ...)-1);
})

setMethodS3("getY", "PlatformDesign", function(this, ...) {
  as.integer(getFeatureInfo(this, "Y", ...)-1);
})

setMethodS3("getOligoSequences", "PlatformDesign", function(this, ...) {
  getFeatureInfo(this, "sequence", ...);
})

setMethodS3("getAlleleAB", "PlatformDesign", function(this, ...) {
  getFeatureInfo(this, "alleleAB", ...);
})

setMethodS3("getUnitNames", "PlatformDesign", function(this, ...) {
  getFeatureInfo(this, "feature_set_name", ...);
})

setMethodS3("isPm", "PlatformDesign", function(this, ...) {
  (getFeatureInfo(this, "feature_type", ...) == "PM");
})


setMethodS3("validate", "PlatformDesign", function(this, subset=1/1000, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the units to be validated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);
  # Currently, AFFX- are non-valid /HB 2006-12-11
  snpUnits <- indexOf(cdf, pattern="^SNP");
  if (is.null(subset)) {
    units <- snpUnits;
  } else if (length(subset) == 1 && subset < 1) {
    subset <- seq(from=1, to=length(snpUnits), length=subset*length(snpUnits));
    subset <- as.integer(subset);
    units <- snpUnits[subset];
  } else {
    units <- subset;
  }

  # Get the CDF-to-PD cell index map
  map <- getReadMap(this);

  verbose && enter(verbose, "Validating ", class(this)[1]);
  verbose && print(verbose, pd);
  verbose && cat(verbose, "Number of units: ", length(units));

  verbose && enter(verbose, "Validating cell indices");
  cellsA <- unlist(getCellIndices(cdf, units=units), use.names=FALSE);
  subset <- map[cellsA];
  cellsO <- getCellIndices(this, subset=subset);
  stopifnot(identical(cellsA, cellsO));
  verbose && exit(verbose);

  verbose && enter(verbose, "Validating unit names");
  unitNamesA <- getUnitNames(cdf, units=units);
  unitNamesO <- getUnitNames(this, subset=subset);
  unitNamesO <- as.character(unique(unitNamesO));
  stopifnot(identical(unitNamesA, unitNamesO));
  verbose && exit(verbose);

  verbose && enter(verbose, "Validating isPm()");
  isPmA <- unlist(isPm(cdf, units=units), use.names=FALSE);
  isPmO <- isPm(this, subset=subset);
  stopifnot(identical(isPmA, isPmO));
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(units);
}, private=TRUE)


############################################################################
# HISTORY:
# 2006-12-11
# o Added validate() to validate consistency between PD package and CDF.
# 2006-12-07
# o getFeatureInfo() now cache to file.  The feature info object is the
#   definitely the biggest object and we rarely need all in memory at the
#   same time.
# o Created in order not to have to load the full platform design package
#   into memory.  Strategy is to split large objects and cache to file.
############################################################################
