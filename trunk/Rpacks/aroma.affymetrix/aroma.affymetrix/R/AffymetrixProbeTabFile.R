###########################################################################/**
# @RdocClass AffymetrixProbeTabFile
#
# @title "The AffymetrixProbeTabFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixProbeTabFile represents an interface to query the data
#  contained in an Affymetrix probe tab file, e.g. 
#  \code{Mapping250K_Nsp_probe_tab}.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of @see "AffymetrixFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# \section{About probe-tab files}{
#  Probe-tab files are provided by Affymetrix and contains information 
#  on the probes.  Note that not necessarily all probes are represented
#  in the files, e.g. typically only PM probes are given and MM are 
#  left to be inferred from the PMs.
#
#  The below is an extract of the \code{Mapping250K_Nsp_probe_tab} file
#  obtained from the Affymetrix website.  Note that columns are separated
#  by TABs.
#  \preformatted{
#   SNP_A-1780270	1633	2398	3	TTGTTAAGCAAGTGAGTTATTTTAT	f	PM	C
#   SNP_A-1780270	1633	2399	3	TTGTTAAGCAAGTGACTTATTTTAT	f	PM	G
#   SNP_A-1780270	1951	1780	-4	GGATAAAATAAAATAACTCACTTGC	r	PM	C
#   ...
#   SNP_A-4241299	2553	1658	4	AAACACATTTTTGGGTCGTAAGGAA	f	PM	G
#  }
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setConstructorS3("AffymetrixProbeTabFile", function(...) {
  extend(AffymetrixFile(..., mustExist=TRUE), "AffymetrixProbeTabFile",
    ".cdf" = NULL,
    "cached:.indexToRowMap" = NULL
  )
}, private=TRUE)


setMethodS3("clearCache", "AffymetrixProbeTabFile", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".indexToRowMap")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)


setMethodS3("as.character", "AffymetrixProbeTabFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  tags <- getTags(this);
  if (!is.null(tags)) {
    s <- paste(s, " Tags: ", paste(tags, collapse=","), ".", sep="");
  }
  s <- c(s, sprintf("Pathname: %s", getPathname(this)));
  s <- c(s, sprintf("File size: %.2fMB", getFileSize(this)/1024^2));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  cdf <- getCdf(this);
  s <- c(s, as.character(cdf));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)

setMethodS3("getChipType", "AffymetrixProbeTabFile", function(this, ...) {
  getChipType(getCdf(this));
})

setMethodS3("getCdf", "AffymetrixProbeTabFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    pattern <- sprintf("_probe_tab$");
    chipType <- gsub(pattern, "", getName(this));
    pathname <- AffymetrixCdfFile$findByChipType(chipType);
    if (is.null(pathname)) {
      throw("Could not located CDF file for chip type pattern: ", pattern);
    }
    cdf <- AffymetrixCdfFile$fromFile(pathname);
    this$.cdf <- cdf;
  }
  cdf;
})


setMethodS3("findByChipType", "AffymetrixProbeTabFile", function(static, chipType, paths=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Searching for probe sequence file");
  verbose && cat(verbose, "Chip type: ", chipType);

  pattern <- paste("_probe_tab$", sep="");
  verbose && cat(verbose, "Filename pattern: ", pattern);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- findAnnotationDataByChipType(chipType, pattern);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # As a backup, search "old" style (code by Ken Simpson)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(pathname)) {
    # Default paths
    paths <- paste(".",
                   getOption("AFFX_SEQUENCE_PATH"),
                   Sys.getenv("AFFX_SEQUENCE_PATH"),
                   "sequence/", "data/sequence/",
                   getOption("AFFX_CDF_PATH"),
                   Sys.getenv("AFFX_CDF_PATH"),
                   "cdf/", "data/cdf/",
                   sep=";", collapse=";");
    pathname <- findFiles(pattern, paths=paths, recursive=TRUE);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # As a backup, search using "old" style v2
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(pathname)) {
    paths <- "annotations";
  
    # First to an exact search
    pattern <- sprintf("^%s_probe_tab$", chipType);
    pathname <- findFiles(pattern=pattern, paths=paths, recursive=TRUE);
    if (length(pathname) == 0) {
      pathname <- NULL;

      # Since Affymetrix is not naming their probe tab files consistently,
      # it might be that they chop of the end of the chip type string.
  
      # 1. Find all probe tab files
      pattern <- sprintf("_probe_tab$");
      pathnames <- findFiles(pattern=pattern, paths=paths, 
                                              firstOnly=FALSE, recursive=TRUE);

      # Found any files?
      if (length(pathnames) > 0) {
        # 2. Extract the part of the filenames containing chip type names
        names <- gsub(pattern, "", basename(pathnames));
      
        # 3. Create patterns out of these
        patterns <- paste("^", names, sep="");
    
        # 4. Keep only those files that match the prefix of the chip type
        keep <- sapply(patterns, function(pattern) {
          (regexpr(pattern, chipType) != -1);
        });
        pathnames <- pathnames[keep];
    
        # 4b. If more than one match was found, keep the longest one
        if (length(pathnames) > 1) {
          names <- names[keep];
          keep <- which.max(nchar(names));
          pathnames <- pathnames[keep];
        }
      }
  
      pathname <- pathnames[1];
    }
  }

  verbose && cat(verbose, "Pathname: ", pathname);

  verbose && exit(verbose);

  pathname;
}, static=TRUE, protected=TRUE)


setMethodS3("fromCdf", "AffymetrixProbeTabFile", function(static, cdf, ...) {
  res <- byChipType(static, getChipType(cdf));
  res$.cdf <- cdf;
  res;
}, static=TRUE);


setMethodS3("fromChipType", "AffymetrixProbeTabFile", function(static, ...) {
  byChipType(static, ...);
}, static=TRUE) 

setMethodS3("byChipType", "AffymetrixProbeTabFile", function(static, chipType, ...) {
  pathname <- AffymetrixProbeTabFile$findByChipType(chipType);
  if (length(pathname) == 0)
    throw("Failed to located the Affymetrix probe tab file: ", chipType);
  newInstance(static, pathname, ...);
}, static=TRUE)


setMethodS3("getIndexToRowMap", "AffymetrixProbeTabFile", function(this, ..., force=FALSE) {
  map <- this$.indexToRowMap;
  if (force || is.null(map)) {
    # c("factor", "integer", "integer", "integer", "character", "factor", "factor", "factor");
    # Read only (X,Y) columns
    colClasses <- rep("NULL", 8);
    colClasses[2:3] <- "integer";
    names(colClasses) <- c("probeSetId", "x", "y", "offset", "sequence", "strand", "type", "allele");
    pathname <- getPathname(this);
    df <- readTable(pathname, colClasses=colClasses, sep="\t", header=FALSE, col.names=names(colClasses));
    
    # Calculate cell indices from (x,y)
    cdf <- getCdf(this);
    indices <- nbrOfColumns(cdf) * df$y + df$x + 1;

    # Get the cell index to (x,y) map.
    map <- rep(NA, nbrOfCells(cdf));
    map[indices] <- seq(along=indices);

    this$.indexToRowMap <- map;
  }

  map;
}, private=TRUE)



setMethodS3("getData", "AffymetrixProbeTabFile", function(this, ...) {
  readDataFrame(this, ...);
})

setMethodS3("readDataFrame", "AffymetrixProbeTabFile", function(this, cells=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Read data");
  verbose && cat(verbose, "Chip type: ", getChipType(this));

  verbose && enter(verbose, "Reading (cell,row) map");
  map <- getIndexToRowMap(this, verbose=less(verbose, 5));
  verbose && cat(verbose, "(cell,row) map:");
  verbose && str(verbose, map);
  verbose && exit(verbose);

  if (is.null(cells)) {
    rows <- map;
  } else {
    rows <- map[cells];
  }

  verbose && cat(verbose, "Unfiltered rows to read:");
  verbose && str(verbose, rows);

  ok <- !is.na(rows);

  if (length(rows[ok]) > 0) {
#    colClasses <- c("factor", "integer", "integer", "integer", "character", "factor", "factor", "factor");
    colClasses <- c("character", "integer", "integer", "integer", "character", "character", "character", "character");
    colClasses[2:3] <- "NULL";
    names(colClasses) <- c("probeSetId", "x", "y", "offset", "sequence", "strand", "type", "allele");
  
    pathname <- getPathname(this);
    verbose && enter(verbose, "Reading data table");
    verbose && cat(verbose, "Pathname: ", pathname);
    verbose && cat(verbose, "Column classes:");
    verbose && str(verbose, as.list(colClasses));
    df <- readTable(pathname, colClasses=colClasses, sep="\t", header=FALSE, col.names=names(colClasses), rows=rows[ok]);
    verbose && exit(verbose);
  } else {
    verbose && cat(verbose, "No rows to read.");
    df <- NULL;
  }

  verbose && enter(verbose, "Expanding result data frame");
  verbose && str(verbose, df);
#  df <- as.list(df);
#  nas <- which(!ok);
#  rm(ok);
#  for (kk in seq(along=df)) {
#    df[[kk]] <- insert(df[[kk]], at=nas, values=NA);
#  }
#  df <- as.data.frame(df);
  colClasses <- sapply(df, FUN=data.class);
  data <- dataFrame(colClasses=colClasses, nrow=length(rows));
  data[ok,] <- df;
  data[!ok,] <- NA;
  rm(ok);
  verbose && str(verbose, data);
  verbose && exit(verbose);

  verbose && exit(verbose);

  data;
}, private=TRUE)


############################################################################
# HISTORY:
# 2008-04-14
# o BUG FIX: readDataFrame() for AffymetrixProbeTabFile would not return
#   the correct number of rows if there were missing cells, which there are.
# o Added verbose output to readDataFrame().
# o Renamed getData() to readDataFrame().
# 2007-06-07
# o When there are no annotation files, findByChipType() of 
#   AffymetrixProbeTabFile would throw "Error in basename(path) : a 
#   character vector argument expected".  Not it returns NULL.
# 2006-12-06
# o Can read probe-tab data from file.  However, for probes not in the file
#   NAs are returned, which typically means that data only for PMs will be
#   available.  It is on the todo list to infer the MM data.  This will
#   require CDF information to match what pair of probes are PM and MM etc.
# o Created.
############################################################################
