###########################################################################/**
# @set "class=ChipEffectSet"
# @RdocMethod extractDataFrame
#
# @title "Extract data as a data.frame for a set of arrays"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @seemethod "extractMatrix".}
#   \item{addNames}{If @TRUE, the first two columns contain the 
#     unit names and the group names according the the CDF, otherwise
#     those two columns are not included.}
#   \item{addUgcMap}{If @TRUE, the columns following the unit and
#     group names contains the (unit, group, cell) index map.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a Jx(2+3+K) @data.frame where J is the number of units, 
#  and K is the number of arrays.  
#  The first two columns, if \code{addNames=TRUE}, contains the 
#  unit names and the group names.
#  The next three columns contains the (unit, group, cell) index map.
#  The last K columns named by the arrays contain the data for the K arrays.
#  No names are set for the rows.
#  The rows are ordered according to \code{units} arguments.
# }
#
# @author
#
# \seealso{
#   @seemethod "extractMatrix".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractDataFrame", "ChipEffectSet", function(this, addNames=FALSE, addUgcMap=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting data for the array set");
  data <- extractMatrix(this, ..., returnUgcMap=TRUE, 
                                                 verbose=less(verbose, 1));

  ugcMap <- attr(data, "unitGroupCellMap");
  attr(data, "unitGroupCellMap") <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  if (addUgcMap) {
    verbose && enter(verbose, "Merging UGC map and extracted data");
    ugcMap <- as.data.frame(ugcMap);
    data <- cbind(ugcMap, data);

    if (addNames) {
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
    }

    verbose && exit(verbose);
  }

  if (addNames) {
    verbose && enter(verbose, "Appending unit and group names from CDF");
    cdf <- getCdf(this);
    verbose && cat(verbose, "CDF chip type: ", 
                                        getChipType(cdf, fullname=TRUE));
    ugNames <- getUnitGroupNamesFromUgcMap(cdf, ugcMap=ugcMap, 
                                              verbose=less(verbose, 10));
    rm(cdf, ugcMap);
    verbose && cat(verbose, "(unit, group) names: ");
    verbose && str(verbose, ugNames);

    ugNames <- as.data.frame(ugNames);
    data <- cbind(ugNames, data);
    rm(ugNames);

    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  data;
}) # extractDataFrame()


############################################################################
# HISTORY:
# 2008-02-12
# o Added arguments 'addNames' and 'addUgcMap'.
# 2008-02-05
# o Created from ChipEffectSet.extractMatrix.R.
############################################################################
