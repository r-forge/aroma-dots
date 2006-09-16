setConstructorS3("GenomeInformation", function(...) {
  this <- extend(AffymetrixFile(...), "GenomeInformation",
    "cached:.data"=NULL
  );
  if (!is.null(getPathname(this)))
    verify(this);
  this;
})

setMethodS3("verify", "GenomeInformation", function(this, ...) TRUE);


setMethodS3("getChipType", "GenomeInformation", function(this, ...) {
  pathname <- getPathname(this);
  basename(dirname(pathname));
})

setMethodS3("fromCdf", "GenomeInformation", function(static, cdf, ...) {
  fromChipType(static, getChipType(cdf), ...);
}, static=TRUE)


setMethodS3("fromChipType", "GenomeInformation", abstract=TRUE);


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
#   as the are ordered in the CDF.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getData", "GenomeInformation", function(this, units=NULL, fields=c("chromosome", "physicalPosition"), orderBy=fields, ...) {
  data <- this$.data;
  if (is.null(data)) {
    # Read the unit names from the corresponding CDF filex
    chipType <- getChipType(this);
    cdfFile <- findCdf(chipType);    
    targetUnitNames <- readCdfUnitNames(cdfFile);

    # Now read the genome information data
    data <- readHg17(this);

    # Reorder genome information data according to CDF file
    idxs <- match(targetUnitNames, data[,1]);
    data <- data[idxs,-1];
    this$.data <- data;
  }

  args <- list(...);
  if (length(args) > 0) {
    for (key in names(args)) {
      value <- args[[key]];
      value <- na.omit(value);
      if (is.function(value)) {
        keep <- value(data[,key]);
      } else {
        keep <- (data[,key] == value);
      }
      data <- data[keep,,drop=FALSE];
    }
    rm(keep);
  }

  if (!is.null(orderBy)) {
    o <- do.call("order", args=as.list(data[,orderBy]));
    data <- data[o,,drop=FALSE];
    rm(o);
  }

  data <- data[,fields,drop=FALSE];
  data;
})


setMethodS3("getUnitIndices", "GenomeInformation", function(this, ..., orderBy=c("chromosome", "physicalPosition"), na.rm=TRUE) {
  df <- getData(this, fields="chromosome", orderBy=orderBy, ...);
  df <- rownames(df);
  df <- as.integer(df);
  if (na.rm)
    df <- df[!is.na(df)];
  df;
})


############################################################################
# HISTORY:
# 2006-09-15
# o Created from DChip.R in old(!) aroma.snp.
# 2005-11-15
# o Added support for skipping header in readSampleInformationFile().
# 2005-10-31
# o Created.
############################################################################  
 