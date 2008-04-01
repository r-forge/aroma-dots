###########################################################################/**
# @set "class=CbsModel"
# @RdocMethod fitOne
#
# @title "Fits the CBS model for one chromosome in one sample"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{data}{}
#   \item{...}{}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the @see "DNAcopy::DNAcopy" object returned by 
#  @see "DNAcopy::segment".
# }
#
# @author
#
# \seealso{
#   Internally @see "DNAcopy::segment" is used.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fitOne", "CbsModel", function(this, data, chromosome, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Fitting CBS");

  args <- list(...);
  keep <- (names(args) %in% names(formals(DNAcopy::segment)));
  fitArgs <- args[keep];

  verbose && enter(verbose, "Setting up CBS data structure");
  nbrOfUnits <- nrow(data);
  chipTypes <- getChipTypes(this);
  data <- DNAcopy::CNA(
    genomdat=data[,"M"], 
    chrom=rep(chromosome, nbrOfUnits),
    data.type="logratio",
    maploc=data[,"x"]
  );
  verbose && str(verbose, data);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calling segment()");
  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && cat(verbose, "Chip types: ", paste(chipTypes, collapse=", "));
  verbose && cat(verbose, "Total number of units: ", nbrOfUnits);
  args <- c(list(data), fitArgs, list(verbose=as.logical(verbose)));
  rm(data, fitArgs);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # segment() writes to stdout; capture it and send it to the verbose object.
  stdout <- capture.output({
    fit <- do.call("segment", args);  # Hmmm... How to do DNAcopy::segment().
  })
  stdout <- paste(stdout, collapse="\n");
  verbose && cat(verbose, stdout);

  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;  
}, private=TRUE) # fitOne()


############################################################################
# HISTORY:
# 2008-03-29
# o fitOne() had argument 'data=data'. Replaced with 'data'.
# 2007-08-20
# o Created from GladModel.fitOne.R.
############################################################################
