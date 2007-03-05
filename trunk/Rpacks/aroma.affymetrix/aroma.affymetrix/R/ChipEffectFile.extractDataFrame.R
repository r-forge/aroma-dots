###########################################################################/**
# @set "class=ChipEffectFile"
# @RdocMethod extractDataFrame
#
# @title "Extract data as a data frame"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{(The subset of units to be matched.
#     If @NULL, all units are considered.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an Jx4 @data.frame where J is the number of units.
#  The names of the last column is the name of the array.
#  The names of the rows are the unit indices (as indexed by the CDF).
#  The rows are ordered according to \code{units} arguments.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractDataFrame", "ChipEffectFile", function(this, ..., units=NULL, field=c("theta", "sdTheta"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  cdf <- getCdf(this);
  if (is.null(units)) {
    nunits <- nbrOfUnits(cdf);
  } else {
    units <- Arguments$getIndices(units, range=c(1,nbrOfUnits(cdf)));
    nunits <- length(units);
  }

  # Argument 'field':
  field <- match.arg(field);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting data for the array");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get cell map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Get unit-to-cell map");
  map <- getCellMap(this, units=units, verbose=less(verbose));
str(map);
  map <- subset(map, ...);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the samples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving data");
  df <- getDataFlat(this, units=map, fields=field, verbose=less(verbose));
  rm(map);
  verbose && exit(verbose);

  verbose && exit(verbose);

  df;
}) # extractDataFrame()

############################################################################
# HISTORY:
# 2007-03-04
# o Created.
############################################################################
