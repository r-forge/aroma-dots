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
  pattern <- sprintf("^%s genome info %s[.]txt$", chipType, version);
  pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE);
  nfiles <- length(pathnames);
  if (nfiles == 0)
    throw("Could not find dChip genome information: ", chipType);
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
    df <- readHg17(this, nrow=10);
  }, error = function(ex) {
    throw("File format error of the dChip genome information file: ", 
                                                  getPathname(this));
  })
  invisible(TRUE);
}, protected=TRUE)


setMethodS3("readHg17", "DChipGenomeInformation", function(this, ..., include=NULL, exclude=c("Expr1002", "Allele A", "dbSNP RS ID")) {
  colClasses <- c(
    "Probe Set ID"="character", 
    "Chromosome"="character",	
    "Expr1002"="integer",	
    "Physical Position"="integer",
    "Allele A"="character",
    "dbSNP RS ID"="character"
  );

  exclude <- setdiff(exclude, include);
  colClasses[names(colClasses) %in% exclude] <- "NULL";

  pathname <- getPathname(this);
  df <- readTable(pathname, ..., colClasses=colClasses, header=TRUE, sep="\t");

  colnames(df) <- toCamelCase(colnames(df));

  df;
}, protected=TRUE)




############################################################################
# HISTORY:
# 2006-09-15
# o Created from DChip.R in old(!) aroma.snp.
# 2005-11-15
# o Added support for skipping header in readSampleInformationFile().
# 2005-10-31
# o Created.
############################################################################  
 