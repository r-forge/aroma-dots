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
#   \item{field}{The field to be extracted.}
#   \item{returnUgcMap}{If @TRUE, the (unit, group, cell) map is returned
#     as an attribute.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an JxK @double @matrix where J is the number of units, 
#  and K is the number of arrays.
#  The names of the columns are the names of the arrays.
#  No names are set for the rows.
#  The rows are ordered according to \code{units} arguments.
# }
#
# @author
#
# \seealso{
#   @seemethod "extractDataFrame".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractMatrix", "ChipEffectSet", function(this, units=NULL, ..., field=c("theta", "sdTheta"), returnUgcMap=FALSE, verbose=FALSE) {
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


  # Settings
  settings <- getOption("aroma.affymetrix.settings");
  gcArrayFrequency <- settings$memory$gcArrayFrequency; 
  if (is.null(gcArrayFrequency))
    gcArrayFrequency <- 10;


  verbose && enter(verbose, "Getting data for the array set");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get cell map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting unit-to-cell map");
  cf <- getFile(this, 1);
  ugcMap <- getUnitGroupCellMap(cf, units=units, verbose=less(verbose));
  ugcMap <- subset(ugcMap, ...);
  verbose && exit(verbose);

  if (nrow(ugcMap) == 0)
    throw("Nothing to return.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate return array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  arrayNames <- getNames(this);
  nbrOfArrays <- length(arrayNames);
  df <- matrix(NA, nrow=nrow(ugcMap), ncol=nbrOfArrays);
  colnames(df) <- arrayNames;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the samples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving sample thetas");
  for (aa in seq_len(nbrOfArrays)) {
    verbose && printf(verbose, "Array %d,\n", aa);
    cf <- getFile(this, aa);
    df[,aa] <- getDataFlat(cf, units=ugcMap, fields=field, 
                                            verbose=less(verbose))[,field];
    if (aa %% gcArrayFrequency == 0) {
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
    }
  } # for (aa in ...)
  verbose && exit(verbose);

  verbose && exit(verbose);

  if (returnUgcMap)
    attr(df, "unitGroupCellMap") <- ugcMap;

  df;
}) # extractMatrix()


############################################################################
# HISTORY:
# 2007-05-26
# o Updated the Rdocs.
# 2007-03-04
# o Created.
############################################################################
