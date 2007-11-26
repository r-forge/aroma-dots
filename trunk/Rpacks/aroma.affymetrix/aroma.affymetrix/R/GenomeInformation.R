###########################################################################/**
# @RdocClass GenomeInformation
#
# @title "The GenomeInformation class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "GenericDataFile".}
#   \item{.verify}{For internal use only.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# @author
#*/###########################################################################
setConstructorS3("GenomeInformation", function(..., .verify=TRUE) {
  extend(GenericDataFile(...), "GenomeInformation",
    "cached:.data"=NULL
  );
})


setMethodS3("as.character", "GenomeInformation", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", this, ...);
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("clearCache", "GenomeInformation", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".data", ".chromosomeStats")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)



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
setMethodS3("verify", "GenomeInformation", function(this, ...) {
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
setMethodS3("getChipType", "GenomeInformation", function(this, ...) {
  # Infer chip type from the first parent directory that has the same name
  # as the chip type of an existing CDF file.
  pathname <- getPathname(this);
  lastPath <- pathname;
  while (TRUE) {
    path <- dirname(lastPath);
    if (path == lastPath)
      break;
    chipType <- basename(path);
    dummy <- AffymetrixCdfFile$findByChipType(chipType);
    if (!is.null(dummy))
      return(chipType)
    lastPath <- path;
  }
  throw("Failed to infer the chip type from the pathname of the genome information file: ", pathname);
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
#   Returns a @see "GenomeInformation" object.
# }
#
# @author
#
# \seealso{
#   @seemethod "fromChipType".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromCdf", "GenomeInformation", function(static, cdf, ...) {
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
#   Returns a @see "GenomeInformation" object.
# }
#
# @author
#
# \seealso{
#   @seemethod "fromCdf".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromChipType", "GenomeInformation", static=TRUE, abstract=TRUE);


setMethodS3("fromDataSet", "GenomeInformation", function(static, dataSet, ...) {
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
#   @seemethod "getUnitIndices".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getData", "GenomeInformation", function(this, units=NULL, fields=c("chromosome", "physicalPosition"), orderBy=NULL, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  data <- this$.data;
  if (is.null(data) || force) {
    verbose && enter(verbose, "Retrieving genome information from file");

    # Read the unit names from the corresponding CDF file
    verbose && enter(verbose, "Reading unit names from CDF file");
    chipType <- getChipType(this);
    cdfFile <- AffymetrixCdfFile$findByChipType(chipType);
    if (is.null(cdfFile))
      throw("Could not located CDF file: ", chipType);
    targetUnitNames <- readCdfUnitNames(cdfFile);
    verbose && exit(verbose);

    # Now read the genome information data
    verbose && enter(verbose, "Reading genome information data");
    data <- readData(this, verbose=less(verbose));
    verbose && str(verbose, data);
    verbose && exit(verbose);

    verbose && enter(verbose, "Reordering units according to the CDF file");
    # Reorder genome information data according to CDF file
    idxs <- match(targetUnitNames, data[,1]);
    data <- data[idxs,,drop=FALSE];
#    data <- data[,-1];
    rownames(data) <- 1:nrow(data);
    verbose && exit(verbose);

    verbose && enter(verbose, "Optimizing default return order");
    # Default ordering
    args <- as.list(data[,fields,drop=FALSE]);
    o <- do.call("order", args=args);
    data <- data[o,,drop=FALSE];
    rm(o);
    verbose && str(verbose, data);
    verbose && exit(verbose);

    if ("chromosome" %in% fields) {
      verbose && enter(verbose, "Replacing 'X' and 'Y' with 23 and 24");
      chr <- data[,"chromosome"];
      chr[chr == "X"] <- 23;
      chr[chr == "Y"] <- 24;
      chr[chr == "Z"] <- 25;
      suppressWarnings({
        chr <- as.integer(chr);
      })
      data[,"chromosome"] <- chr;
      rm(chr);
      verbose && str(verbose, data);
      verbose && exit(verbose);
    }

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
    rr <- match(units, rownames(data));
    data <- data[rr,,drop=FALSE];
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


setMethodS3("getUnitsOnChromosome", "GenomeInformation", function(this, chromosomes, region=NULL, ...) {
  allChromosomes <- getChromosomes(this);

  # Argument 'chromosomes':
  chromosomes <- Arguments$getChromosomes(chromosomes, range=range(allChromosomes), ...);

  # Argument 'region':
  if (!is.null(region)) {
    region <- Arguments$getNumerics(region);
    if (length(region) != 2) {
      throw("Argument 'region' must be a numeric vector of length 2: ", 
                                                         length(region));
    }
    if (region[1] > region[2]) {
      throw("Argument 'region' is not ordered: c(", 
                                      paste(region, collapse=", "), ")");
    }
  }

  data <- getData(this);
  keep <- (data[,"chromosome"] %in% chromosomes);
  data <- data[keep,,drop=FALSE];

  # Stratify by physical position?
  if (!is.null(region)) {
    pos <- data[,"physicalPosition"];
    keep <- (!is.na(pos) & (region[1] <= pos & pos <= region[2]));
    data <- data[keep,,drop=FALSE];
  }

  # Extract the units
  units <- rownames(data);
  units <- as.integer(units);

  units;
})


#setMethodS3("readData", "GenomeInformation", abstract=TRUE);

setMethodS3("readData", "GenomeInformation", function(this, ...) {
  readTableInternal(this, ...);
}, protected=TRUE);


setMethodS3("readTableInternal", "GenomeInformation", function(this, pathname, colClasses=NULL, ..., include=NULL, exclude=NULL, fill=TRUE, verbose=FALSE) {
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
  
  df <- readTable(pathname, colClasses=colClasses, header=TRUE, sep="\t", fill=fill, ..., verbose=less(verbose));

  colnames(df) <- toCamelCase(colnames(df));
  verbose && exit(verbose);

  df;
}, private=TRUE)




setMethodS3("nbrOfUnits", "GenomeInformation", function(this, ...) {
  data <- getData(this, fields=1);
  nrow(data);
})

setMethodS3("getDataColumns", "GenomeInformation", function(this, ...) {
  data <- getData(this);
  colnames(data);
}, private=TRUE)


###########################################################################/**
# @RdocMethod getUnitIndices
#
# @title "Gets unit indices ordered along the chromosome"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getData".}
#  \item{na.rm}{If @TRUE, non-defined unit indices are excluded, otherwise
#      not.}
# }
#
# \value{
#   Returns an @integer @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "getData".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getUnitIndices", "GenomeInformation", function(this, ..., na.rm=TRUE) {
  df <- getData(this, fields="chromosome", ...);
  df <- rownames(df);
  suppressWarnings({  # Suppress warnings about NAs
    df <- as.integer(df);
  })
  if (na.rm)
    df <- df[!is.na(df)];
  df;
})


###########################################################################/**
# @RdocMethod getPositions
#
# @title "Gets the physical positions for a set of units"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getData".}
#  \item{na.rm}{If @TRUE, non-defined unit indices are excluded, otherwise
#      not.}
# }
#
# \value{
#   Returns an @integer @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "getData".
#   @seemethod "getUnitIndices".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPositions", "GenomeInformation", function(this, ..., na.rm=FALSE) {
  df <- getData(this, fields="physicalPosition", ...);
  suppressWarnings({  # Suppress warnings about NAs
    df <- as.integer(df[,1]);
  })
  if (na.rm)
    df <- df[!is.na(df)];
  df;
})



setMethodS3("getChromosomes", "GenomeInformation", function(this, ..., force=FALSE) {
  chromosomes <- this$.chromosomes;
  if (is.null(chromosomes) || force) {
    chromosomes <- unique(getData(this, fields="chromosome")[,1]);
    chromosomes <- sort(chromosomes);

##    # Sort in order of (1:22,"X","Y")
##    chromosomeMap <- c(1:22,"X","Y", NA);
##    o <- match(chromosomeMap, chromosomes);
##    chromosomes <- chromosomes[o];
##    chromosomes <- chromosomes[!is.na(chromosomes)];
    
    this$.chromosomes <- chromosomes;
  }
  chromosomes;
})


setMethodS3("getChromosomeStats", "GenomeInformation", function(this, na.rm=TRUE, ..., force=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Retrieving chromosome statistics");

  stats <- this$.chromosomeStats;
  if (is.null(stats) || force) {
    chromosomes <- getChromosomes(this);
    nbrOfChromosomes <- length(chromosomes);
    stats <- matrix(NA, nrow=nbrOfChromosomes, ncol=3);
    colnames(stats) <- c("min", "max", "n");
    rownames(stats) <- chromosomes;
    for (chr in chromosomes) {
      verbose && enter(verbose, sprintf("Chromosome %d of %d", 
                                               chr, length(chromosomes)));
      pos <- getPositions(this, chromosome=chr);
      r <- range(pos, na.rm=na.rm);
      stats[chr,1:2] <- r;
      stats[chr,3] <- length(pos);
      verbose && exit(verbose);
    }
    this$.chromosomeStats <- stats;
  } else {
    verbose && cat(verbose, "Found cached results");
  }

  verbose && exit(verbose);

  stats;  
})


###########################################################################/**
# @RdocMethod plotDensity
#
# @title "Plots the density of SNPs for a given chromosome"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chromosome}{The chromsome to be displayed.}
#  \item{...}{Additional arguments passed to @seemethod "getPositions".}
#  \item{adjust}{A bandwidth parameter for the density estimator.}
#  \item{xlab}{The label on the x-axis.  If @NULL, a default will generated.}
#  \item{main}{The title of the plot.  If @NULL, a default will generated.}
#  \item{annotate}{If @TRUE, the plot is annotated with extra information.}
# }
#
# \value{
#   Returns (invisibly) the estimated density.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("plotDensity", "GenomeInformation", function(this, chromosome, ..., adjust=1/20, xlab=NULL, main=NULL, annotate=TRUE) {
  # Get the positions of these units
  pos <- getPositions(this, chromosome=chromosome, ...);

  # Estimate the SNP density
  d <- density(pos/10^6, from=0, adjust=adjust);

  # Plot the density
  if (is.null(xlab))
    xlab <- "Physical position (Mb)";
  if (is.null(main))
    main <- sprintf("SNP density on chromsome %s", chromosome);
  plot(d, lwd=2, xlab=xlab, main=main);

  # Annotate?
  if (annotate) {
    stext(sprintf("%d SNPs", d$n), side=3, pos=1, line=-1, cex=0.8);
    cdf <- AffymetrixCdfFile$fromChipType(getChipType(this));
    stextChipType(cdf);
  }

  invisible(d);
})


############################################################################
# HISTORY:
# 2007-11-25
# o Now getUnitsOnChromosome() of GenomeInformation can indentify units from
#   multiple chromosomes.
# 2007-10-30
# o Now argument 'chromosome' for getUnitsOnChromosome() needs to be 
#   specified explictly. Before its default was '23'.
# 2007-09-10
# o Now readData() in GenomeInformation is no longer abstract, but tries
#   naively to read data using readTableInternal() as is.  That will make
#   it slightly easier to add new genome information files.
# 2007-08-12
# o BUG GIX: clearCache() would not clear the genome stats cache.
# o BUG FIX: Subsetting with getData(..., chromosome=...) would return NAs
#   for units with missing information.
# 2007-06-11
# o BUG FIX: Used non-existing 'gi' instead of 'this' in plotDensity() of
#   GenomeInformation.
# 2007-03-15
# o Updated GenomeInformation to return chromosomes as indices and never
#   with 'X' and 'Y' regardless of source.  This is part of a moving the
#   package to handle chromosomes as indices so that it will work with
#   any genome.
# 2007-02-28
# o Added argument 'region' to getUnitsOnChromosome().
# 2007-01-25
# o BUG FIX: Added 'fill=TRUE' to readTableInternal().
# 2007-01-06
# o Renamed getFields() to getDataColumns().
# 2006-11-29
# o Added getUnitsOnChromosome().
# 2006-09-16
# o Added plotDensity(), getChromosomes(), getChromosomeStats() etc.
# o Improved getData().  Updated getUnitIndices() and getPositions().
# 2006-09-15
# o Added Rdoc comments.
# o Created from DChip.R in old(!) aroma.snp.
# 2005-11-15
# o Added support for skipping header in readSampleInformationFile().
# 2005-10-31
# o Created.
############################################################################  
 
