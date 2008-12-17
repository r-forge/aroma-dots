###########################################################################/**
# @RdocClass HaarSegModel
#
# @title "The HaarSegModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Haar wavelet-based segmentation (HaarSeg) 
#  model [1]. 
# }
# 
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "ChipEffectSetTuple".}
#   \item{breaksFdrQ}{Default tuning parameters specific to the HaarSeg 
#         algorithm.}
#   \item{...}{Arguments passed to the constructor of 
#              @see "CopyNumberSegmentationModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#  @see "CopyNumberSegmentationModel".
# }
#
# \references{
#   [1] Ben-Yaacov E. and Eldar YC. \emph{A fast and flexible method for the segmentation of aCGH data}, Bioinformatics, 2008.
#   \url{http://www.ee.technion.ac.il/Sites/People/YoninaEldar/Info/software/HaarSeg.htm}
# }
#*/########################################################################### 
setConstructorS3("HaarSegModel", function(cesTuple=NULL, breaksFdrQ=0.0001, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(cesTuple)) {
    require("HaarSeg") || throw("Package not loaded: HaarSeg");
  }

  # Argument 'breaksFdrQ':
  breaksFdrQ <- Arguments$getDouble(breaksFdrQ, range=c(0,1));

  extend(CopyNumberSegmentationModel(cesTuple=cesTuple, ...), "HaarSegModel",
    .breaksFdrQ = breaksFdrQ
  )
})


setMethodS3("getAsteriskTags", "HaarSegModel", function(this, collapse=NULL, ...) {
  tags <- "HAAR";

  # Add class-specific tags
  if (isPaired(this))
    tags <- c(tags, "paired");

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)


setMethodS3("fitOne", "HaarSegModel", function(this, data, chromosome, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Fitting HaarSeg");

  args <- list(...);
  keep <- (names(args) %in% names(formals(HaarSeg::haarSeg)));
  fitArgs <- args[keep];

  verbose && enter(verbose, "Setting up HaarSeg data structure");
  nbrOfUnits <- nrow(data);
  chipTypes <- getChipTypes(this);

  # Order data along chromosome
  o <- order(data[,"x"]);
  data <- data[o,,drop=FALSE];
  rm(o);

  # Setup arguments to HaarSeg::haarSeg()
  args <- list(
    I = data[,"M"],
    breaksFdrQ = this$.breaksFdrQ
  );
  verbose && str(verbose, args);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calling HaarSeg::haarSeg()");
  verbose && cat(verbose, "Chip types: ", paste(chipTypes, collapse=", "));
  verbose && cat(verbose, "Total number of units: ", nbrOfUnits);
  verbose && cat(verbose, "Additional arguments:");
  verbose && str(verbose, fitArgs);
  args <- c(args, fitArgs);
  rm(fitArgs);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # segment() writes to stdout; capture it and send it to the verbose object.
  verbose && cat(verbose, "Arguments to segmentation method:");
  verbose && str(verbose, args);
  verbose && enter(verbose, "Calling HaarSeg::haarSeg() function");
  stdout <- capture.output({
    fit <- do.call("haarSeg", args);
  })
  rm(args);
  verbose && str(verbose, fit);
  verbose && exit(verbose);
  stdout <- paste(stdout, collapse="\n");
  verbose && cat(verbose, stdout);

  verbose && exit(verbose);

  verbose && enter(verbose, "Setting up return HaarSeg object");
  fit <- list(
    output=fit, 
    data=list(M=data[,"M"], x=data[,"x"], chromosome=chromosome)
  );
  class(fit) <- "HaarSeg";
  rm(data);
  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;  
}, private=TRUE) # fitOne()




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# HaarSeg class wrappers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("extractRawCopyNumbers", "HaarSeg", function(object, ...) {
  data <- object$data;
  RawCopyNumbers(cn=data$M, x=data$x, chromosome=data$chromosome);
})


setMethodS3("drawCnRegions", "HaarSeg", function(object, ...) {
  output <- object$output;
  regions <- output$SegmentsTable;
  data <- object$data;
  x <- data$x;

  # Identify the locus indices where the regions starts and ends
  starts <- regions[,1];
  counts <- regions[,2];
  ends <- starts+counts-1;

  # Translate to genomic positions
  starts <- x[starts];
  ends <- x[ends];

  # Get the mean levels of each region
  means <- regions[,3];

  nbrOfSegments <- length(means);
  xx <- cbind(starts, ends);
  yy <- cbind(means, means);
  for (rr in 1:nrow(xx)) {
    lines(x=xx[rr,], y=yy[rr,], ...);
  };
})



setMethodS3("extractCopyNumberRegions", "HaarSeg", function(object, ...) {
  output <- object$output;
  regions <- output$SegmentsTable;
  data <- object$data;
  x <- data$x;

  # Identify the locus indices where the regions starts and ends
  starts <- regions[,1];
  counts <- regions[,2];
  ends <- starts+counts-1;

  # Translate to genomic positions
  starts <- x[starts];
  ends <- x[ends];

  # Get the mean levels of each region
  means <- regions[,3];

  CopyNumberRegions(
    chromosome=data$chromosome,
    start=starts, 
    stop=ends, 
    mean=means,
    count=counts
  );
})



##############################################################################
# HISTORY:
# 2007-12-17
# o Now using the HaarSeg package (put together by HB).
# 2007-12-16
# o Created.
############################################################################## 
