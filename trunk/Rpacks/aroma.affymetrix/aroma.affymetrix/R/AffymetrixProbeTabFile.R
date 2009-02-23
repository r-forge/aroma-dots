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
  extend(TabularTextFile(..., mustExist=TRUE),
                  c("AffymetrixProbeTabFile", uses("AromaPlatformInterface")),
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

  s <- NextMethod("as.character", this, ...);
  class <- class(s);

#  s <- sprintf("%s:", class(this)[1]);
#  s <- c(s, sprintf("Name: %s", getName(this)));
#  tags <- getTags(this);
#  if (!is.null(tags)) {
#    s <- paste(s, " Tags: ", paste(tags, collapse=","), ".", sep="");
#  }
#  s <- c(s, sprintf("Pathname: %s", getPathname(this)));
#  s <- c(s, sprintf("File size: %.2fMB", getFileSize(this)/1024^2));
#  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  cdf <- getCdf(this);
  s <- c(s, as.character(cdf));

  class(s) <- class;
  s;
}, private=TRUE)


setMethodS3("translateFullName", "AffymetrixProbeTabFile", function(this, name, ...) {
  name <- gsub("[._]probe(|[._]tab)", "", name);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Turning special Affymetrix tags into regular tags
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  name <- gsub("[.](CN)(|[,_].*)", ",\\1\\2", name);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Patching incorrect Affymetrix file names
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  map <- c(
    "HG-U133-PLUS" = "HG-U133_Plus_2",
    "Mapping50K_Hind" = "Mapping50K_Hind240",
    "Mapping50K_Xba" = "Mapping50K_Xba240"
  );

  patterns <- sprintf("^(%s)(,.*|$)", names(map));
  idx <- sapply(patterns, FUN=regexpr, name);
  idx <- which(idx != -1);
  if (length(idx) > 0) {
    idx <- idx[1];
    name <- gsub(patterns[idx], sprintf("%s\\2", map[idx]), name);
  }
  
  name <- NextMethod("translateFullName", this, name, ...);
  name;
})


setMethodS3("hasColumnHeader", "AffymetrixProbeTabFile", function(this, ...) {
  # Infer if there is a column header or not
  hdr <- readRawHeader(this);
  hasColumnHeader <- any(regexpr("^(PROBESET|Probe Set)", hdr$topRows[[1]]) != -1);
  hasColumnHeader;
}, protected=TRUE)



setMethodS3("translateColumnNames", "AffymetrixProbeTabFile", function(this, names, ...) {
##  # Remove all prefixes "PROBE_"
##  names <- gsub("^PROBE_", "", names);

  # Convert 'Foo_baR_dOo' and 'FOO.baR.dOo' to 'Foo baR dOo'?
  if (any(regexpr("[_.]", names) != -1)) {
    names <- gsub("[_.]", " ", names);

    # Convert to lower case
    names <- tolower(names);
  }


  lcNames <- tolower(names);
  names[lcNames == "probe set name"] <- "probeset id";
  names[lcNames == "probe x"] <- "probe x pos";
  names[lcNames == "probe y"] <- "probe y pos";

  # Finally, convert 'foo bar doo' to 'fooBarDoo'
  names <- toCamelCase(names);

  names <- NextMethod("translateColumnNames", this, names, ...);
  names;
}, protected=TRUE)



setMethodS3("getColumnNames", "AffymetrixProbeTabFile", function(this, ...) {
  columns <- NextMethod("getColumnNames", this, ...);

  if (is.null(columns)) {
    # Has column header?
    topRow <- getHeader(this)$topRows[[1]];
    nbrOfColumns <- length(topRow);
    defColumns <- sprintf("V%02d", 1:nbrOfColumns);
    columns <- rep(NA, nbrOfColumns);

    if (hasColumnHeader(this)) {
      columns <- topRow;
    } else {
      # Identify PROBESET_ID
      pattern <- "^[a-zA-Z]+[a-zA-Z0-9]_[a-zA-Z0-9].*";
      pattern <- toAsciiRegExprPattern(pattern);
      idx <- grep(pattern, topRow);
      if (length(idx) > 0) {
        columns[idx] <- "PROBESET_ID";
      }
  
      # Identify PROBE_SEQUENCE
      pattern <- "^[ACGT]{25}$";
      idx <- grep(pattern, topRow);
      if (length(idx) > 0)
        columns[idx] <- "PROBE_SEQUENCE";
  
      # Identify PROBE_TYPE
      pattern <- "^(PM|MM)$";
      idx <- grep(pattern, topRow);
      if (length(idx) > 0)
        columns[idx] <- "PROBE_TYPE";
  
      # Identify TARGET_STRANDEDNESS
      pattern <- "^(\\+|-|f|r)$";
      idx <- grep(pattern, topRow);
      if (length(idx) > 0)
        columns[idx] <- "TARGET_STRANDEDNESS";
  
      # Identify ALLELE
      pattern <- "^[ACGT]$";
      idx <- grep(pattern, topRow);
      if (length(idx) > 0)
        columns[idx] <- "ALLELE";
  
      # Guess remaining
      idxs <- which(is.na(columns));
      nidxs <- length(idxs);
      if (nidxs >= 1)
        columns[idxs[1]] <- "PROBE_X_POS";
      if (nidxs >= 2)
        columns[idxs[2]] <- "PROBE_Y_POS";
      if (nidxs >= 3)
        columns[idxs[3]] <- "PROBE_INTERROGATION_POSITION";
    }
  }

  # Finally, translate any column names
  columns <- translateColumnNames(this, columns);

  columns;
})


setMethodS3("getChipType", "AffymetrixProbeTabFile", function(this, ...) {
  pattern <- sprintf("[._]probe_tab$");
  chipType <- gsub(pattern, "", getName(this));

  # Patch non-consistent Affymetrix filenames
  if (chipType == "Mapping10K") {
    if (getFileSize(this) != 14452965) {
      chipType <- "Mapping10K_Xba142";
    } else {
      chipType <- "Mapping10K_Xba131";
    }
  }


  chipType;
})


setMethodS3("getCdf", "AffymetrixProbeTabFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    chipType <- getChipType(this);
    pathname <- AffymetrixCdfFile$findByChipType(chipType);
    if (is.null(pathname))
      throw("Could not located CDF file for chip type pattern: ", chipType);
    cdf <- AffymetrixCdfFile$fromFile(pathname);
    this$.cdf <- cdf;
  }
  cdf;
})


setMethodS3("findByChipType", "AffymetrixProbeTabFile", function(static, chipType, what=c("", "CN"), paths=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  if (is.null(what)) {
    what <- "";
  } else {
    what <- match.arg(what);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Searching for probe sequence file");
  verbose && cat(verbose, "Chip type: ", chipType);

  pattern <- paste("[._]probe[._]tab$", sep="");
  if (what == "CN") {
    pattern <- paste("[.]CN", pattern, sep="");
  }
  verbose && cat(verbose, "Filename pattern: ", pattern);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- findAnnotationDataByChipType(chipType, pattern, 
                                                          firstOnly=FALSE);
  if (length(pathname) > 1) {
    verbose && cat(verbose, "Located files:");
    verbose && print(verbose, pathname);

    # Identify the shortest matching filename
    filenames <- basename(pathname);
    pattern2 <- sprintf("^(%s)(.*)(%s)", chipType, pattern);
    tags <- gsub(pattern2, "\\2", filenames);
    idx <- which.min(nchar(tags));
    pathname <- pathname[idx];
  }
  verbose && print(verbose, pathname);

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
    pattern <- sprintf("^%s[._]probe_tab$", chipType);
    pathname <- findFiles(pattern=pattern, paths=paths, recursive=TRUE);
    if (length(pathname) == 0) {
      pathname <- NULL;

      # Since Affymetrix is not naming their probe tab files consistently,
      # it might be that they chop of the end of the chip type string.
  
      # 1. Find all probe tab files
      pattern <- "[._]probe_tab$";
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
  res <- byChipType(static, chipType=getChipType(cdf), nbrOfCells=nbrOfCells(cdf), ...);
  res$.cdf <- cdf;
  res;
}, static=TRUE);


setMethodS3("fromChipType", "AffymetrixProbeTabFile", function(static, ...) {
  byChipType(static, ...);
}, static=TRUE) 

setMethodS3("byChipType", "AffymetrixProbeTabFile", function(static, chipType, what=NULL, nbrOfCells=NULL, ...) {

  pathname <- AffymetrixProbeTabFile$findByChipType(chipType, what=what, ...);
  if (length(pathname) == 0)
    throw("Failed to located the Affymetrix probe tab file: ", chipType);
  res <- newInstance(static, pathname, ...);

  # Validate?
  if (!is.null(nbrOfCells)) {
    # Possible? /HB 2009-02-10
  }

  res;
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


setMethodS3("readDataFrame2", "AffymetrixProbeTabFile", function(this, cells=NULL, ..., verbose=FALSE) {
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
# 2008-07-07
# o Updated findByChipType() to search for pattern "[._]probe[._]tab$", 
#   because all combinations exist.
# 2008-07-06
# o Added argument 'what' to findByChipType() and byChipType().
# o Updated translateColumnNames() to handle more file format versions.
# o Updated translateFullName() to fix more inconsistent filenames.
# 2008-06-30
# o Added getChipType() patch for infering the chiptype from inconsistent
#   Affymetrix filenames.
# 2008-06-12
# o Added translateFullName() to fix incorrect Affymetrix filenames, e.g.
#   'Mapping50K_Hind_probe_tab' instead of 'Mapping50K_Hind240_probe_tab'.
# o Added default filename pattern to ".*[._]probe_tab$".
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
