###########################################################################/**
# @RdocClass UgpGenomeInformation
#
# @title "The UgpGenomeInformation class"
#
# \description{
#  @classhierarchy
#
#  This class represents Aroma UGP genome information files.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "GenomeInformation".}
#   \item{.ugp}{For internal use only.}
#   \item{.verify}{For internal use only.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# @author
#*/###########################################################################
setConstructorS3("UgpGenomeInformation", function(..., .ugp=NULL, .verify=TRUE) {
  this <- extend(GenomeInformation(..., .verify=FALSE), "UgpGenomeInformation",
    .ugp = .ugp
  );
  if (.verify) {
    if (!is.null(getPathname(this)))
      verify(this);
  }
  this;
})


setMethodS3("getAromaUgpFile", "UgpGenomeInformation", function(this, ..., force=FALSE) {
  ugp <- this$.ugp;
  if (force || is.null(ugp)) {
    ugp <- AromaUgpFile(getPathname(this), ...);
    this$.ugp <- ugp;
  }
  ugp;
}, protected=TRUE);


setMethodS3("findByChipType", "UgpGenomeInformation", function(static, ...) {
  AromaUgpFile$findByChipType(...);
}, static=TRUE, protected=TRUE)


###########################################################################/**
# @RdocMethod fromChipType
#
# @title "Defines a UgpGenomeInformation object by chip type"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chipType}{A @character string.}
#  \item{version}{An optional @character string specifying the version
#    string, if more than one version is available.}
#  \item{...}{Not used.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "UgpGenomeInformation" object.  
#  If no file was not found, an error is thrown.
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
setMethodS3("fromChipType", "UgpGenomeInformation", function(static, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  ugp <- AromaUgpFile$fromChipType(...);
  pathname <- getPathname(ugp);

  verbose && enter(verbose, "Instantiating ", class(static)[1]);
  verbose && cat(verbose, "Pathname: ", pathname);

  res <- newInstance(static, filename=pathname, path=NULL, .ugp=ugp);
  verbose && print(verbose, res);
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
})


setMethodS3("verify", "UgpGenomeInformation", function(this, ...) {
  tryCatch({
    df <- readData(this, nrow=10);
  }, error = function(ex) {
    throw("File format error of the UGP genome information file (", 
                                 ex$message, "): ", getPathname(this));
  })
  invisible(TRUE);
}, private=TRUE)


setMethodS3("readData", "UgpGenomeInformation", function(this, nrow=NULL, ..., verbose=FALSE) {
  verbose && enter(verbose, "Reading data from UGP file");

  ugp <- getAromaUgpFile(this);
  verbose && print(verbose, ugp, level=-20);

  if (is.null(nrow)) {
    verbose && cat(verbose, "Reading all ", nbrOfUnits(ugp), " units");
    res <- ugp[,,drop=FALSE];
  } else {
    units <- 1:nrow;
    verbose && cat(verbose, "Reading ", length(units), " units");
    res <- ugp[units,,drop=FALSE];
  }

  colnames(res) <- c("chromosome", "physicalPosition");
  verbose && exit(verbose);

  res;
})

setMethodS3("getDataColumns", "UgpGenomeInformation", function(this, ...) {
  c("chromosome", "physicalPosition");
}, private=TRUE)


setMethodS3("getData", "UgpGenomeInformation", function(this, units=NULL, fields=getDataColumns(this), orderBy=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  data <- this$.data;
  if (is.null(data) || force) {
    verbose && enter(verbose, "Retrieving genome information from file");

    # Now read the genome information data
    ugp <- getAromaUgpFile(this);
    cc <- match(fields, getDataColumns(this));
    missing <- fields[is.na(cc)];
    if (length(missing)) {
      throw("Unknown fields: ", paste(missing, collapse=", "));
    }
  
    verbose && enter(verbose, "Reading genome information data");
    data <- ugp[,,drop=FALSE];
    colnames(data) <- getDataColumns(this);
    verbose && str(verbose, data);
    verbose && exit(verbose);

    # Store in cache
    this$.data <- data;

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    
    verbose && exit(verbose);
  }

  # Subset by unit?
  if (!is.null(units)) {
    # Map the unit indicies to the row names
    data <- data[units,,drop=FALSE];
  }

  # Stratify by field values?
  args <- list(...);
  if (length(args) > 0) {
    for (key in names(args)) {
      # Get the values to be stratified upon.
      values <- data[,key,drop=FALSE];

      # Get the test (value or function)
      test <- args[[key]];
      test <- na.omit(test);
      if (is.function(test)) {
        keep <- test(values);
      } else {
        keep <- (values == test);
        keep <- (keep & !is.na(keep));
      }
      data <- data[keep,,drop=FALSE];
    }
    rm(keep);
  }

  # Reorder?
  if (!is.null(orderBy)) {
    o <- do.call("order", args=as.list(data[,orderBy,drop=FALSE]));
    data <- data[o,,drop=FALSE];
    rm(o);
  }

  # Extract a subset of fields?
  if (!is.null(fields))
    data <- data[,fields, drop=FALSE];

  data;
})


setMethodS3("getChromosomes", "UgpGenomeInformation", function(this, force=FALSE, ...) {
  chromosomes <- this$.chromosomes;

  if (is.null(chromosomes) || force) {
    ugp <- getAromaUgpFile(this);
    chromosomes <- ugp[,1];
    gc <- gc();
    chromosomes <- unique(chromosomes);
    gc <- gc();
    chromosomes[chromosomes == 0] <- NA;
    chromosomes <- chromosomes[!is.na(chromosomes)];
    chromosomes <- sort(chromosomes);
    chromosomes <- as.integer(chromosomes);
    this$.chromosomes <- chromosomes;
  }

  chromosomes;
})


############################################################################
# HISTORY:
# 2007-09-11
# o Created.
############################################################################  
