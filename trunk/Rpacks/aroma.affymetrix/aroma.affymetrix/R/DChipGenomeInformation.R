###########################################################################/**
# @RdocClass DChipGenomeInformation
#
# @title "The DChipGenomeInformation class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "GenomeInformation".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# @author
#*/###########################################################################
setConstructorS3("DChipGenomeInformation", function(...) {
  this <- extend(GenomeInformation(...), "DChipGenomeInformation");
  if (!is.null(getPathname(this)))
    verify(this);
  this;
})

setMethodS3("fromChipType", "DChipGenomeInformation", function(static, chipType, path="annotations", version=NULL, ...) {
  # Argument 'path' & 'chipType':
  path <- filePath(path, chipType, expandLinks="any");
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Search for genome information files
  if (is.null(version))
    version <- ".*";

  # Create a filename pattern
  pattern <- sprintf("^.*( |_)genome( |_)info(| |_).*%s[.](txt|xls)$", version);

  pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE);
  nfiles <- length(pathnames);
  if (nfiles == 0) {
    throw("Could not find dChip genome information: ", chipType);
  }

  # More than one match?
  if (nfiles > 1) {
    pathnames <- sort(pathnames);
    warning("Found more than one matching dChip genome information file, but returning only the last one: ", paste(pathnames, collapse=", "));
    pathnames <- rev(pathnames);
  }
  pathname <- pathnames[1];
   
  newInstance(static, pathname);
})


setMethodS3("verify", "DChipGenomeInformation", function(this, ...) {
  tryCatch({
    df <- readData(this, nrow=10);
  }, error = function(ex) {
    throw("File format error of the dChip genome information file: ", 
                                                  getPathname(this));
  })
  invisible(TRUE);
}, private=TRUE)


setMethodS3("readData", "DChipGenomeInformation", function(this, ...) {
  readFcns <- list(
    "^Mapping10K" = read50KHg17,
    "^Mapping50K" = read50KHg17,
    "^Mapping250K" = read250KHg17
  );

  chipType <- getChipType(this);

  # Try to read with the designated read function.
  res <- NULL;
  for (kk in seq(along=readFcns)) {
    pattern <- names(readFcns)[kk];
    if (regexpr(pattern, chipType) != -1) {
      readFcn <- readFcns[[kk]];
      tryCatch({
        res <- readFcn(this, ...);
      }, error=function(ex) {})
    }
  }

  # If failed, re-try using all read functions.
  if (is.null(res)) {
    for (kk in seq(along=readFcns)) {
      readFcn <- readFcns[[kk]];
      tryCatch({
        res <- readFcn(this, ...);
      }, error=function(ex) {})
      if (!is.null(res))
        break;
    }
  }

  if (is.null(res)) {
    throw("Cannot read dChip annotation data.  No predefined read function available for this chip type: ", chipType);
  }

  res;
})




setMethodS3("read250KHg17", "DChipGenomeInformation", function(this, ..., exclude=c("Expr1002", "Allele A", "dbSNP RS ID")) {
  colClasses <- c(
    "Probe Set ID"="character", 
    "Chromosome"="character",	
    "Expr1002"="integer",	
    "Physical Position"="integer",
    "Allele A"="character",
    "dbSNP RS ID"="character"
  );
  readTableInternal(this, pathname=getPathname(this), colClasses=colClasses, exclude=exclude, ...);
}, private=TRUE)


setMethodS3("read50KHg17", "DChipGenomeInformation", function(this, ..., exclude=c("Genetic Map", "Strand", "Allele A", "dbSNP RS ID", "Freq AfAm", "heterrate")) {
  colClasses <- c(
    "Probe Set ID"="character", 
    "Chromosome"="character",	
    "Physical Position"="integer",
    "Genetic Map"="",               # Often empty?!
    "Strand"="character",
    "Allele A"="character",
    "dbSNP RS ID"="character",
    "Freq AfAm"="double",
    "heterrate"="double"            # Often empty?!
  );

  readTableInternal(this, pathname=getPathname(this), colClasses=colClasses, exclude=exclude, ...);
}, private=TRUE)



############################################################################
# HISTORY:
# 2007-01-22
# o Made fromChipType() identify genome information files more robustly, 
#   e.g. the exact chip type does not have to be part of the prefix.
#   Indeed, it now accepts the default dChip filenames for the 10K, the
#   50K and the 250K SNP chips.
# o Made readData() to support also unknown chip types. That is, if a chip
#   type is not among the hardwired ones, the method will still try to read
#   it using one of the known read functions.  For instance, the genome info
#   file for 10K chips have the same format as the one for the 100K chips.
# 2006-09-15
# o Created from DChip.R in old(!) aroma.snp.
# 2005-11-15
# o Added support for skipping header in readSampleInformationFile().
# 2005-10-31
# o Created.
############################################################################  
