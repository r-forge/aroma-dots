###########################################################################/**
# @RdocClass DChipSnpInformation
#
# @title "The DChipSnpInformation class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "SnpInformation".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# @author
#*/###########################################################################
setConstructorS3("DChipSnpInformation", function(...) {
  this <- extend(SnpInformation(...), "DChipSnpInformation");
  if (!is.null(getPathname(this)))
    verify(this);
  this;
})

setMethodS3("fromChipType", "DChipSnpInformation", function(static, chipType, path="annotations", version=NULL, ...) {
  # Argument 'path' & 'chipType':
  path <- filePath(path, chipType, expandLinks="any");
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Search for SNP information files
  if (is.null(version))
    version <- ".*";

  # Create a filename pattern
  pattern <- sprintf("^%s snp info[.]txt$", chipType, version);

  pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE);
  nfiles <- length(pathnames);
  if (nfiles == 0) {
    throw("Could not find dChip SNP information: ", chipType);
  }

  if (nfiles > 1) {
    pathnames <- sort(pathnames);
    warning("Found more than one matching dChip SNP information file, but returning only the last one: ", paste(pathnames, collapse=", "));
    pathnames <- rev(pathnames);
  }
  pathname <- pathnames[1];
   
  newInstance(static, pathname);
})


setMethodS3("verify", "DChipSnpInformation", function(this, ...) {
  tryCatch({
    df <- readData(this, nrow=10);
  }, error = function(ex) {
    throw("File format error of the dChip SNP information file: ", 
                                                  getPathname(this));
  })
  invisible(TRUE);
}, protected=TRUE)


setMethodS3("readData", "DChipSnpInformation", function(this, ...) {
  chipType <- getChipType(this);
  if (regexpr("^Mapping50K", chipType) != -1) {
    read50K(this, ...);
  } else if (regexpr("^Mapping250K", chipType) != -1) {
    read250K(this, ...);
  } else {
    throw("Cannot read dChip SNP information file.  No predefined read function available for this chip type: ", chipType);
  }
})




setMethodS3("read250K", "DChipSnpInformation", function(this, ..., exclude=c("dbSNP RS ID", "Flank", "FreqAsia", "FreqAfAm", "FreqCauc")) {
  # Example with TABs replaced by semicolons:
  # Probe Set ID;dbSNP RS ID;Flank;Fragment Length Start Stop;FreqAsia;FreqAfAm;FreqCauc
  # SNP_A-1780520;rs16994928;ggatagtgttgacctc[A/G]agtacaggtttcaaaa;496 // 47873735 // 47874230;0.0 ;0.11;0.0 
  colClasses <- c(
    "Probe Set ID"="character", 
    "dbSNP RS ID"="character",
    "Flank"="character",	
    "Fragment Length Start Stop"="character",	
    "FreqAsia"="double",
    "FreqAfAm"="double",
    "FreqCauc"="double"
  );
  readTableInternal(this, pathname=getPathname(this), colClasses=colClasses, exclude=exclude, ...);
}, protected=TRUE)

setMethodS3("read50K", "DChipSnpInformation", function(this, ..., exclude=c("dbSNP RS ID", "Flank", "FreqAsian", "FreqAfAm", "FreqCauc")) {
  colClasses <- c(
    "Probe Set ID"="character", 
    "dbSNP RS ID"="character",
    "Flank"="character",	
    "Fragment Length Start Stop"="character",	
    "FreqAsian"="double",
    "FreqAfAm"="double",
    "FreqCauc"="double"
  );
  readTableInternal(this, pathname=getPathname(this), colClasses=colClasses, exclude=exclude, ...);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2006-09-17
# o Created from DChipGenomeInformation.R.
############################################################################  
