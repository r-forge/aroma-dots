###########################################################################/**
# @set "class=ChipEffectSet"
# @RdocMethod extractMatrix
#
# @title "Extract data as a matrix for a set of arrays"
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
#  Returns an Jx(3+K) @double @matrix where J is the number of units, 
#  and K is the number of arrays.
#  The names of the columns (after the 3rd) are the names of the arrays.
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
setMethodS3("extractMatrix", "ChipEffectSet", function(this, ..., units=NULL, field=c("theta", "sdTheta"), verbose=FALSE) {
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


  verbose && enter(verbose, "Getting data for the array set");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get cell map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Get unit-to-cell map");
  cf <- getFile(this, 1);
  map <- getCellMap(cf, units=units, verbose=less(verbose));
  map <- subset(map, ...);
  verbose && exit(verbose);

  if (nrow(map) == 0)
    throw("Nothing to return.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate return array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  arrayNames <- getNames(this);
  nbrOfArrays <- length(arrayNames);
  df <- matrix(NA, nrow=nrow(map), ncol=3+nbrOfArrays);
  colnames(df) <- c(colnames(map), arrayNames);
  df[,1:3] <- as.matrix(map);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the samples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving sample thetas");
  for (aa in seq_len(nbrOfArrays)) {
    cf <- getFile(this, aa);
    values <- getDataFlat(cf, units=map, fields=field, verbose=less(verbose))[,field];
    df[,3+aa] <- values;
    rm(values);
  } # for (aa in ...)
  verbose && exit(verbose);

  verbose && exit(verbose);

  df;
}) # extractMatrix()


############################################################################
# HISTORY:
# 2007-03-04
# o Created.
############################################################################
