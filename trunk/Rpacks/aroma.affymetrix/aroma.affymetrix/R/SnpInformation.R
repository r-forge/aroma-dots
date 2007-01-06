###########################################################################/**
# @RdocClass SnpInformation
#
# @title "The SnpInformation class"
#
# \description{
#  @classhierarchy
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
# @author
#*/###########################################################################
setConstructorS3("SnpInformation", function(...) {
  this <- extend(AffymetrixFile(...), "SnpInformation",
    "cached:.data"=NULL
  );
  if (!is.null(getPathname(this)))
    verify(this);
  this;
})


###########################################################################/**
# @RdocMethod verify
#
# @title "Verifies the correctness of the underlying file"
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
#   Returns (visibly) @TRUE if the file is valid, otherwise an error is
#   thrown.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("verify", "SnpInformation", function(this, ...) {
  TRUE;
}, private=TRUE)



###########################################################################/**
# @RdocMethod getChipType
#
# @title "Gets the chip type of this genome information set"
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
#   Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChipType", "SnpInformation", function(this, ...) {
  pathname <- getPathname(this);
  basename(dirname(pathname));
})


###########################################################################/**
# @RdocMethod fromCdf
#
# @title "Static method to define a genome information set from a CDF"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{cdf}{An @see "AffymetrixCdfFile".}
#  \item{...}{Additional argument passed to @seemethod "fromChipType".}
# }
#
# \value{
#   Returns a @see "SnpInformation" object.
# }
#
# @author
#
# \seealso{
#   @seemethod "fromChipType".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromCdf", "SnpInformation", function(static, cdf, ...) {
  fromChipType(static, getChipType(cdf), ...);
}, static=TRUE)




###########################################################################/**
# @RdocMethod fromChipType
#
# @title "Static method to define a genome information set by chip type"
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
#   Returns a @see "SnpInformation" object.
# }
#
# @author
#
# \seealso{
#   @seemethod "fromCdf".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromChipType", "SnpInformation", abstract=TRUE);


setMethodS3("fromDataSet", "SnpInformation", function(static, dataSet, ...) {
  cdf <- getCdf(dataSet);
  chipType <- getChipType(cdf);
  fromChipType(static, chipType, ...);
}, static=TRUE)



###########################################################################/**
# @RdocMethod getData
#
# @title "Gets all or a subset of the genome information data"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{units}{The units for which the data should be returned.}
#  \item{fields}{The fields to be returned.}
#  \item{orderBy}{The fields by which the returned data frame should be 
#      ordered.}
#  \item{...}{Named arguments used to select a subset of the units to be
#      returned.  Either a value to be compared to or a @function returning
#      @TRUE or @FALSE.}
# }
#
# \value{
#   Returns a @data.frame, where the row names correspond to unit indices
#   as defined by the CDF.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getData", "SnpInformation", function(this, units=NULL, fields=c("fragmentLength", "start", "stop"), orderBy=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  data <- this$.data;
  if (is.null(data) || force) {
    verbose && enter(verbose, "Requiring SNP information from file");

    # Read the unit names from the corresponding CDF file
    verbose && enter(verbose, "Reading unit names from CDF file");
    chipType <- getChipType(this);
    cdfFile <- findCdf(chipType);
    if (is.null(cdfFile))
      throw("Could not located CDF file: ", chipType);
    targetUnitNames <- readCdfUnitNames(cdfFile);
    verbose && exit(verbose);

    # Now read the SNP information data
    verbose && enter(verbose, "Reading SNP information data");
    data <- readData(this, verbose=less(verbose));
    verbose && exit(verbose);

    verbose && enter(verbose, "Splitting up the fragment length, start & stop details");
    cc <- which("fragmentLengthStartStop" == colnames(data));
    lss <- data[,cc,drop=TRUE];
    nas <- (lss == "---");
    lss[nas] <- "-//-//-";
    lss <- strsplit(lss, split="//");
    len <- sapply(lss, FUN=length);
    if (any(len != 3))
      throw("Internal error: Unrecognized FLSS.");
    lss <- unlist(lss, use.names=FALSE);
    lss <- as.integer(lss);
    lss <- matrix(lss, ncol=3, byrow=TRUE);
    rownames(lss) <- 1:nrow(lss);
    colnames(lss) <- c("fragmentLength", "start", "stop");
    data <- cbind(data[,-cc,drop=FALSE], lss);
    rm(lss);
    verbose && exit(verbose);

    verbose && enter(verbose, "Reordering units according to the CDF file");
    idxs <- match(targetUnitNames, data[,1]);
    data <- data[idxs,];
    rownames(data) <- 1:nrow(data);
#    data <- data[idxs,-1];
    verbose && exit(verbose);
    # Store in cache
    this$.data <- data;

    verbose && exit(verbose);
  }

  # Subset by unit?
  if (!is.null(units)) {
    # Map the unit indicies to the row names
    rr <- match(units, rownames(data));
    data <- data[rr,,drop=TRUE];
  }
  # Stratify by field values?
  args <- list(...);
  if (length(args) > 0) {
    for (key in names(args)) {
      # Get the values to be stratified upon.
      values <- data[,key];

      # Get the test (value or function)
      test <- args[[key]];
      test <- na.omit(test);
      if (is.function(test)) {
        keep <- test(values);
      } else {
        keep <- (values == test);
      }
      data <- data[keep,,drop=FALSE];
    }
    rm(keep);
  }

  # Reorder?
  if (!is.null(orderBy)) {
    o <- do.call("order", args=as.list(data[,orderBy]));
    data <- data[o,,drop=FALSE];
    rm(o);
  }

  # Extract a subset of fields?
  if (!is.null(fields))
    data <- data[,fields, drop=FALSE];

  data;
})


setMethodS3("nbrOfUnits", "SnpInformation", function(this, ...) {
  data <- getData(this, fields=1);
  nrow(data);
})

setMethodS3("getFields", "SnpInformation", function(this, ...) {
  data <- getData(this);
  colnames(data);
})


setMethodS3("readData", "SnpInformation", abstract=TRUE);

setMethodS3("readTableInternal", "SnpInformation", function(this, pathname, colClasses=NULL, ..., include=NULL, exclude=NULL, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Reading tabular data from file");
  pathname <- getPathname(this);
  verbose && cat(verbose, "Pathname: ", pathname);

  verbose && cat(verbose, "Argument 'include': ", paste(include, collapse=", "));
  verbose && cat(verbose, "Argument 'exclude': ", paste(exclude, collapse=", "));
  exclude <- setdiff(exclude, include);
  verbose && cat(verbose, "exclude\\include: ", paste(exclude, collapse=", "));
  colClasses[names(colClasses) %in% exclude] <- "NULL";
  toRead <- names(colClasses)[colClasses != "NULL"];
  verbose && cat(verbose, "Columns to be read: ", paste(toRead, collapse=", "));
  
  df <- readTable(pathname, colClasses=colClasses, header=TRUE, sep="\t", ..., verbose=less(verbose));

  colnames(df) <- toCamelCase(colnames(df));
  verbose && exit(verbose);

  df;
}, private=TRUE)


setMethodS3("getFragmentLengths", "SnpInformation", function(this, ...) {
  data <- getData(this, ..., fields="fragmentLength");
  fl <- data[,1];
  fl <- as.integer(fl);
  fl;
})


setMethodS3("getFragmentStarts", "SnpInformation", function(this, ...) {
  data <- getData(this, ..., fields="start");
  fl <- data[,1];
  fl <- as.integer(fl);
  fl;
})


setMethodS3("getFragmentStops", "SnpInformation", function(this, ...) {
  data <- getData(this, ..., fields="stop");
  fl <- data[,1];
  fl <- as.integer(fl);
  fl;
})



############################################################################
# HISTORY:
# 2006-09-17
# o Created from GenomeInformation.R.
############################################################################  
