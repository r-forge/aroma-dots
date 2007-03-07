###########################################################################/**
# @RdocClass GenotypeCallFile
#
# @title "The GenotypeCallFile class"
#
# \description{
#  @classhierarchy
#
#  The abstract GenotypeCallFile class represents a file containing genotype
#  calls for a given genotype parameter estimates.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{File format}{
#  A GenotypeCallFile object is stored to file as a
#  @see "R.huge::FileByteVector".
# }
# 
# @author
#
# @visibility "private"
#*/###########################################################################
setConstructorS3("GenotypeCallFile", function(..., cdf=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(cdf)) {
    if (!inherits(cdf, "AffymetrixCdfFile"))
      throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  this <- extend(AffymetrixFile(...), "GenotypeCallFile",
    .cdf = cdf
  );

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})


setMethodS3("as.character", "GenotypeCallFile", function(this, ...) {
  s <- NextMethod("as.character");
  s <- c(s, "CDF:");
  cdf <- getCdf(this);
  s <- c(s, as.character(cdf));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)



###########################################################################/**
# @RdocMethod create
#
# @title "Static method creating a genotype call file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file.}
#   \item{path}{An optional path to the file.}
#   \item{cdf}{An @see "AffymetrixCdfFile".}
#   \item{overwrite}{If @TRUE, an already existing file is overwritten,
#      otherwise not.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "GenotypeCallFile".
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
setMethodS3("create", "GenotypeCallFile", function(static, filename, path=NULL, cdf, overwrite=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' & 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path, 
                                                    mustNotExist=!overwrite);

  if (!inherits(cdf, "AffymetrixCdfFile")) {
    throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Creates the genotype file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  header <- sprintf("chipType=%s;", getChipType(cdf));
  comments <- sprintf("header=%s\n", header);
  calls <- FileByteVector(pathname, length=nbrOfUnits(cdf), comments=comments);
  on.exit(close(calls));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Defines the genotype file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  newInstance(static, pathname, cdf=cdf, ...);
}, static=TRUE, private=TRUE)



###########################################################################/**
# @RdocMethod fromFile
#
# @title "Static method defining a genotype call file object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file.}
#   \item{path}{An optional path to the file.}
#   \item{...}{Not used.}
#   \item{.checkArgs}{If @TRUE, arguments are validated and parsed, 
#     otherwise not.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "GenotypeCallFile".
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
setMethodS3("fromFile", "GenotypeCallFile", function(static, filename, path=NULL, ..., .checkArgs=TRUE, verbose=FALSE) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    # Argument 'filename' & 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, 
                                                         mustExist=TRUE);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }
  } else {
    pathname <- filename;
  }

  # Try to define class
  res <- newInstance(static, pathname, ...);

  # Try to get file header
  header <- getHeader(res);

  # Make sure there is a chip type attribute
  chipType <- header$chipType;
  if (is.null(chipType)) {
    throw("File format error: Genotype call file does not contain header attribute 'chipType': ", pathname);
  }

  res;
}, static=TRUE)


setMethodS3("getHeader", "GenotypeCallFile", function(this, ...) {
  pathname <- getPathname(this);
  fdata <- FileByteVector(pathname);
  on.exit(close(fdata));
  comments <- getComments(fdata);
  header <- gsub("header=(.*)\n.*", "\\1", comments);
  header <- unlist(strsplit(header, split=";"));
  header <- strsplit(header, split="=");
  names <- sapply(header, FUN=.subset, 1);
  header <- sapply(header, FUN=.subset, 2);
  names(header) <- names;
  header <- as.list(header);
  header;
}, private=TRUE)



###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF"
#
# \description{
#  @get "title" for this call file.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCdfFile".
# }
#
# \details{
#  The CDF is looked up from the chip type attribute part of the file
#  header.
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
setMethodS3("getCdf", "GenotypeCallFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    header <- getHeader(this);
    chipType <- header$chipType;
    cdf <- AffymetrixCdfFile$fromChipType(chipType);
    this$.cdf <- cdf;
  }
  cdf;
})

setMethodS3("setCdf", "GenotypeCallFile", function(this, cdf, ...) {
  # Argument 'cdf':
  if (!inherits(cdf, "AffymetrixCdfFile"))
    throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
 
  this$.cdf <- cdf;
})


# setMethodS3("getChipType", "GenotypeCallFile", function(this, ...) {
#   cdf <- getCdf(this);
#   chipType <- getChipType(cdf, fullname=FALSE);
#   chipType;
# })



###########################################################################/**
# @RdocMethod "["
#
# @title "Gets the genotype calls"
#
# \description{
#  @get "title" for a subset of units.
# }
#
# @synopsis
#
# \arguments{
#   \item{i}{An @integer @vector specifying the subset of units to query.
#     If @NULL, all units are considered.}
#   \item{drop}{If @TRUE, single dimensions are dropped.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @factor @vector with levels \code{-}, \code{AA}, \code{AB},
#  \code{BB}, and \code{NC}.
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
setMethodS3("[", "GenotypeCallFile", function(this, i=NULL, drop=TRUE, ...) {
  res <- readUnits(this, units=i, ...);
  if (drop && length(res) == 1)
    res <- drop(res);
  res;
})


###########################################################################/**
# @RdocMethod "[<-"
#
# @title "Sets the genotype calls"
#
# \description{
#  @get "title" for a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{i}{An @integer @vector specifying the subset of units to update.
#     If @NULL, all units are considered.}
#   \item{value}{A @vector of the same length as \code{i} with genotype
#     calls.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   Internally @seemethod "updateUnits" is used.
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("[<-", "GenotypeCallFile", function(this, i=NULL, value) {
  updateUnits(this, units=i, calls=value);
})



###########################################################################/**
# @RdocMethod readUnits
#
# @title "Gets the genotype calls"
#
# \description{
#  @get "title" for a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{An @vector specifying the subset of units to update.
#     If @NULL, all units are considered.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @factor @vector with levels \code{-}, \code{AA}, \code{AB},
#  \code{BB}, and \code{NC}.
# }
#
# @author
#
# \seealso{
#   @seemethod "updateUnits".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("readUnits", "GenotypeCallFile", function(this, units=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  cdf <- getCdf(this);
  units <- convertUnits(cdf, units=units, keepNULL=TRUE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Open data file
  pathname <- getPathname(this);
  fdata <- FileByteVector(pathname);
  on.exit(close(fdata));

  # Retrieve data
  calls <- fdata[units];

  # Turn into a factor
  labels <- c("-", "AA", "AB", "BB", "NC");
  calls <- factor(calls, levels=0:4, labels=labels);

  calls;
})



###########################################################################/**
# @RdocMethod updateUnits
#
# @title "Updates the genotype calls"
#
# \description{
#  @get "title" for a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{An @vector specifying the subset of units to update.
#     If @NULL, all units are considered.}
#   \item{calls}{A @vector of the same length as \code{units} with 
#     genotype calls.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "readUnits".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("updateUnits", "GenotypeCallFile", function(this, units=NULL, calls, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  cdf <- getCdf(this);
  units <- convertUnits(cdf, units=units);

  # Argument 'calls':
  if (is.numeric(calls)) {
    calls <- as.integer(calls);
  } else if (is.factor(calls)) {
    calls <- as.integer(calls);
  } else if (is.character(calls)) {
    calls <- match(calls, c("AA", "AB", "BB", "NC"), nomatch=0);
  }

  if (length(calls) != length(units)) {
    if (length(calls) == 1) {
      calls <- rep(calls, length.out=length(units));
    } else {
      throw("Argument 'units' and 'calls' are of different lengths: ", 
                                     length(calls), " != ", length(units));
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Open data file
  pathname <- getPathname(this);
  fdata <- FileByteVector(pathname);
  on.exit(close(fdata));

  # Update data
  fdata[units] <- calls;

  invisible(this);
})



###########################################################################/**
# @RdocMethod importFromCrlmmFile
#
# @title "Imports genotype calls from a CRLMM tab-delimited call file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file.}
#   \item{path}{An optional path to the file.}
#   \item{sample}{The name of the sample column to be extracted/imported.}
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns (invisibly) the number of units imported.
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
setMethodS3("importFromCrlmmFile", "GenotypeCallFile", function(this, filename="calls.txt", path=NULL, sample=getName(this), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing genotype calls from CRLMM output file.");
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Parse external genotype call file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Pathname: ", pathname);

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
  header <- gsub("[.](c|C)(e|E)(l|L)$", "", header);
  nbrOfSamples <- length(header);

  verbose && printf(verbose, "Detected %d samples.", nbrOfSamples);

  if (is.character(sample)) {
    sampleIdx <- grep(sample, header);
    if (length(sampleIdx) == 0) {
      throw("Could not import genotype calls. No sample '", sample, 
                                          "' in call file: ", pathname);
    } else if (length(sampleIdx) > 2) {
      throw("Could not import genotype calls. More than one match for sample '", sample, "' in call file: ", pathname);
    }
    sample <- sampleIdx;
  }
  verbose && printf(verbose, "Requested sample in column %d.", sample);
  sample <- Arguments$getIndex(sample, range=c(1, nbrOfSamples));

  verbose && enter(verbose, "Reading data");
  colClasses <- c("character", rep("NULL", nbrOfSamples));
  colClasses[sample+1] <- "integer";
  calls <- read.table(fh, colClasses=colClasses, row.names=1, sep=sep, header=FALSE);
  ncalls <- nrow(calls);
  verbose && cat(verbose, "Number of rows read: ", ncalls);
  verbose && exit(verbose);

  # Sort calls according to CDF
  verbose && enter(verbose, "Sort according to CDF file");
  cdf <- getCdf(this);
  verbose && cat(verbose, "Number of units in CDF: ", nbrOfUnits(cdf));
  units <- indexOf(cdf, names=rownames(calls));
  nbad <- which(is.na(units));
  
  # Check for units not in CDF and exclude them
  if (length(nbad) == 0) {
    verbose && cat(verbose, "All ", length(units), " SNPs read are available in the CDF file.");
  } else {
    warning(sprintf("%d (%%.1f) of the SNPs in the call file was not found in the CDF file: %s", nbad, 100*nbad/ncalls, pathname));
    units <- units[-nbad];
  }
  verbose && exit(verbose);

  calls <- calls[,1];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update call file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  updateUnits(this, units=units, calls=calls, verbose=less(verbose));

  verbose && exit(verbose);

  invisible(length(units));
}, private=TRUE)



############################################################################
# HISTORY:
# 2007-01-15
# o Added setCdf().
# o Added more Rdoc comments.
# 2006-12-13
# o Revived.
# 2006-10-01
# o Can now import a single CRLMM column.
# 2006-09-30
# o Created.
############################################################################
