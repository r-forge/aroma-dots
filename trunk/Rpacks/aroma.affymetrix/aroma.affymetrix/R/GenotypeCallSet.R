###########################################################################/**
# @RdocClass GenotypeCallSet
#
# @title "The GenotypeCallSet class"
#
# \description{
#  @classhierarchy
#
#  The GenotypeCallSet class represents a set of genotype-call files.
# }
# 
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "GenotypeCallFile":s.}
#   \item{...}{Arguments passed to @see "AffymetrixFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
# @visibility "private"
#*/###########################################################################
setConstructorS3("GenotypeCallSet", function(files=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    lapply(files, FUN=function(df) {
      if (!inherits(df, "GenotypeCallFile"))
        throw("Argument 'files' contains a non-GenotypeCallFile object: ", class(df));
    })
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }

  extend(AffymetrixFileSet(files=files, ...), "GenotypeCallSet"
  )
})


###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the genotype call set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("as.character", "GenotypeCallSet", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Number of arrays: %d", nbrOfArrays(this)));
  s <- c(s, sprintf("Total file size: %.2fMB", getFileSize(this)/1024^2));
  s <- c(s, "CDF:");
  cdf <- getCdf(this);
  s <- c(s, as.character(cdf));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)



###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF structure for this CEL set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCdfFile" object.
# }
#
# @author
#
# \seealso{
#   @seemethod "setCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCdf", "GenotypeCallSet", function(this, ...) {
  getCdf(this$files[[1]], ...);
})


# setMethodS3("getChipType", "GenotypeCallSet", function(this, ...) {
#   cdf <- getCdf(this);
#   chipType <- getChipType(cdf, fullname=FALSE);
#   chipType;
# })


###########################################################################/**
# @RdocMethod setCdf
#
# @title "Sets the CDF structure for this genotype call set"
#
# \description{
#  @get "title".  This structures is assigned to all genotype call files
#    in the set.
# }
#
# @synopsis
#
# \arguments{
#   \item{cdf}{An @see "AffymetrixCdfFile" object.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "getCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("setCdf", "GenotypeCallSet", function(this, cdf, ...) {
  # Nothing to do?
  oldCdf <- getCdf(this);
  if (equals(cdf, oldCdf))
    return(invisible(this));

  # Set the CDF for all CEL files
  lapply(this, setCdf, cdf, ...);

  # Have to clear the cache 
  clearCache(this);

  invisible(this);
})



###########################################################################/**
# @RdocMethod fromFiles
#
# @title "Static method defining a set of genotype call files"
#
# \description{
#  @get "title" in a directory.
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The directory where the call files are.}
#   \item{patttern}{The filename pattern used to identify call files.}
#   \item{...}{Arguments passed to \code{fromFiles()} of 
#      @see "AffymetrixFileSet".}
#   \item{fileClass}{The @see "R.oo::Class" of each genotype file.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "GenotypeCallSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("fromFiles", "GenotypeCallSet", function(static, path, pattern="[.]calls$", ..., fileClass="GenotypeCallFile", verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # S3 method dispatch does not work for static methods.
  ds <- fromFiles.AffymetrixFileSet(static, path=path, pattern=pattern, ..., fileClass=fileClass, verbose=less(verbose));

  # Use the same CDF object for all CEL files.
  setCdf(ds, getCdf(ds));

  ds;
}, static=TRUE)



setMethodS3("getFullName", "GenotypeCallSet", function(this, parent=1, ...) {
  NextMethod("getFullName", this, parent=parent, ...);
})


###########################################################################/**
# @RdocMethod nbrOfArrays
#
# @title "Gets the number of arrays in the file set"
#
# \description{
#   @get "title".
#   This is just a wrapper for \code{nbrOfFiles()}.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfArrays", "GenotypeCallSet", function(this, ...) {
  nbrOfFiles(this, ...);
})


###########################################################################/**
# @RdocMethod as.GenotypeCallSet
# @alias as.GenotypeCallSet.list
# @alias as.GenotypeCallSet.default
#
# @title "Coerce an object to an GenotypeCallSet object"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Other arguments passed to @see "base::list.files".}
# }
#
# \value{
#   Returns an @see "GenotypeCallSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.GenotypeCallSet", "GenotypeCallSet", function(object, ...) {
  object;
})

setMethodS3("as.GenotypeCallSet", "list", function(object, ...) {
  GenotypeCallSet(object, ...);
})

setMethodS3("as.GenotypeCallSet", "default", function(object, ...) {
  throw("Cannot coerce object to an GenotypeCallSet object: ", mode(object));
})



###########################################################################/**
# @RdocMethod "["
# @aliasmethod "[["
#
# @title "Gets the genotype calls"
#
# \description{
#  @get "title" for a subset of units and a subset of arrays.
# }
#
# @synopsis
#
# \arguments{
#   \item{i}{An @integer @vector specifying the subset of units to query.
#     If @NULL, all units are considered.}
#   \item{j}{An @integer @vector specifying the subset of arrays to query.
#     If @NULL, all arrays are considered.}
#   \item{...}{Additional arguments passed to @seemethod "readUnits".}
#   \item{drop}{If @TRUE, single dimensions are dropped.}
# }
#
# \value{
#  Returns a @factor @matrix (@vector) 
#  with levels \code{-}, \code{AA}, \code{AB}, \code{BB}, and \code{NC}.
# }
#
# @author
#
# \seealso{
#   Internally @seemethod "readUnits" is used.
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("[", "GenotypeCallSet", function(this, i=NULL, j=NULL, ..., drop=FALSE) {
  res <- readUnits(this, units=i, arrays=j, ...);
  if (drop && length(res) == 1)
    res <- res[[1]];
  res;
})

setMethodS3("[[", "GenotypeCallSet", function(this, units=NULL, ...) {
  this[units=units, ..., drop=TRUE];
})


setMethodS3("readUnits", "GenotypeCallSet", function(this, units=NULL, arrays=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  cdf <- getCdf(this);
  units <- convertUnits(cdf, units=units, keepNULL=TRUE);
  if (is.null(units)) {
    nbrOfUnits <- nbrOfUnits(cdf);
  } else {
    nbrOfUnits <- length(units);
  }

  # Argument 'arrays':
  if (is.null(arrays)) {
    arrays <- seq(this);
  } else if (is.character(arrays)) {
    arrayNames <- arrays;
    arrays <- match(arrayNames, getNames(this));
    if (any(is.na(arrays))) {
      throw("Argument 'arrays' contains unknown array names: ", 
                      paste(arrayNames[is.na(arrays)], collapse=", "));
    }
  } else {
    arrays <- Arguments$getIndices(arrays, range=c(1, nbrOfFiles(this)));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="readUnits", class=class(this)[1], 
              units=units, arrays=arrays);
  id <- digest(key);
  res <- this$.readUnitsCache[[id]];
  if (!is.null(res)) {
    verbose && cat(verbose, "readUnits(): Returning cached data");
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the pathnames of all files
  pathnames <- getPathnames(this)[arrays];
  nbrOfArrays <- length(arrays);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read data from file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Retrieving genotype calls for ", nbrOfArrays, " arrays");

  calls <- matrix("-", nrow=nbrOfUnits, ncol=nbrOfArrays);

  for (kk in seq(length=nbrOfArrays)) {
    array <- arrays[kk];
    df <- getFile(this, array);
    verbose && enter(verbose, sprintf("Array %s ('%s') ", kk, getName(df)));
    value <- readUnits(df, units=units, ...);
    # If more than one column, assume first is the call column
    if (!is.null(ncol(value)))
      value <- value[,1];
    calls[,kk] <- as.character(value);
    verbose && exit(verbose);
  }
  rm(value, df, array);

  # Turn into a factor
  dim <- dim(calls);  # ...which will drop dimensions
  calls <- as.factor(calls);
  dim(calls) <- dim;
  colnames(calls) <- getNames(this)[arrays];
  rownames(calls) <- units;

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "readUnits(): Updating cache");
  this$.readUnitsCache <- list();
  this$.readUnitsCache[[id]] <- calls;

  calls;
})


setMethodS3("createFromCrlmmFile", "GenotypeCallSet", function(this, filename="calls.txt", path=NULL, outPath=filePath("calls", getChipType(cdf, fullname=FALSE)), cdf, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);

  # Argument 'cdf':
  if (!is.null(cdf)) {
    if (is.character(cdf)) {
      cdf <- AffymetrixCdfFile$fromChipType(cdf);
    }

    if (!inherits(cdf, "AffymetrixCdfFile"))
      throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  # Argument 'outPath':
  outPath <- Arguments$getWritablePath(outPath);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Creating genotype-call set from CRLMM call file.");
  verbose && cat(verbose, "Source pathname: ", pathname);
  verbose && cat(verbose, "Destination path: ", outPath);

  # Read column header
  fh <- file(pathname, open="r");
  on.exit(close(fh));

  # Get the sample names
  header <- readLines(con=fh, n=1);
  sep <- " ";
  if (regexpr("\t", header) != -1)
    sep <- "\t";
  header <- unlist(strsplit(header, split=sep));
  header <- trim(header);
  header <- gsub("^\"", "", header);
  header <- gsub("\"$", "", header);
  samples <- gsub("[.](c|C)(e|E)(l|L)$", "", header);
  nbrOfSamples <- length(samples);
  verbose && printf(verbose, "Detected %d samples: %s", 
                              nbrOfSamples, paste(samples, collapse=", "));

  # Create all sample files
  verbose && enter(verbose, "Creating ", nbrOfSamples, " genotype-call files");
  files <- vector("list", nbrOfSamples);
  for (kk in seq(along=samples)) {
    filename <- sprintf("%s.calls", samples[kk]);
    files[[kk]] <- GenotypeCallFile$create(filename, path=outPath, cdf=cdf);
  }
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading data");
  colClasses <- c("character", rep("integer", nbrOfSamples));
  calls <- read.table(fh, colClasses=colClasses, row.names=1, sep=sep, header=FALSE);
  ncalls <- nrow(calls);
  verbose && cat(verbose, "Number of rows read: ", ncalls);
  verbose && exit(verbose);

  # Sort calls according to CDF
  verbose && enter(verbose, "Sort according to CDF file");
  verbose && cat(verbose, "Number of units in CDF: ", nbrOfUnits(cdf));
  units <- indexOf(cdf, names=rownames(calls));
  nbad <- which(is.na(units));
  
  # Check for units not in CDF and exclude them
  if (length(nbad) == 0) {
    verbose && cat(verbose, "All ", length(units), " SNPs read are available in the CDF file.");
  } else {
    warning(sprintf("%d (%%.1f) of the SNPs in the call file was not found in the CDF file: %s", nbad, 100*nbad/ncalls, pathname));
    units <- units[-nbad];
    calls <- calls[-nbad,];
  }
  verbose && exit(verbose);

  # Update
  verbose && enter(verbose, "Updating all call files");
  for (kk in 1:ncol(calls)) {
    files[[kk]][units] <- calls[,kk];
  }
  verbose && exit(verbose);

  verbose && exit(verbose);

  newInstance(this, files=files);
}, private=TRUE)



############################################################################
# HISTORY:
# 2006-12-17
# o BUG FIX: readUnits() would not accept names in the 'arrays' argument.
# 2006-12-14
# o Updated readUnits() to return a matrix of factors instead.
# 2006-10-01
# o Can now import a single CRLMM column.
# o Created.
############################################################################
