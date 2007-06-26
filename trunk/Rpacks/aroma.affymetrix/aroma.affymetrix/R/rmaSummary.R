###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod rmaSummary
#
# @title "Calculates the RMA expression summary"
#
# \description{
#  @get "title".
#
# Models the observed signal as normal background + exponential signal, then
# applies quantile normalisation followed by a summary by either median
# polish or robust regression.
#
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The path to save the expression summaries.}
#   \item{name}{Name of the data set containing the expression summaries.}
#   \item{bgPath}{Directory in which to store background-adjusted signals.}
#   \item{normPath}{Directory in which to store normalised signals.}
#   \item{summaryMethod}{Either "medianpolish" or "rlm".}
#   \item{overwrite}{If @TRUE, already adjusted arrays are overwritten,
#     unless skipped, otherwise an error is thrown.}
#   \item{skip}{If @TRUE, the array is not normalized if it already exists.}
#   \item{...}{Additional arguments passed to some of the internal functions.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet" containing the expression summaries.
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#
# \seealso{
#  See package \pkg{affy}.
#  @seeclass
# }
#*/###########################################################################
setMethodS3("rmaSummary", "AffymetrixCelSet", function(this, path=NULL, name="rma", bgPath=NULL, normPath=NULL, summaryMethod="rlm", ..., verbose=FALSE) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  if (is.null(path)) {
    # Path structure: /rma/<data set name>/chip_data/<chip type>/
    path <- file.path(name, getName(this), "chip_data", getChipType(cdf));
  }
  if (!is.null(path)) {
    # Verify this path (and create if missing)
    path <- Arguments$getWritablePath(path);
  }

  if (identical(getPath(this), path)) {
    throw("Cannot compute expression summaries. Argument 'path' refers to the same path as the path of the raw probe level data: ", path);
  }
  mkdirs(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # do the background correction step

  verbose && enter(verbose, "Performing background correction");

  dsBG <- bgAdjustRma(this, path=bgPath, name=name, ..., verbose=verbose);

  verbose && exit(verbose);

  verbose && enter(verbose, "Normalising");
  
  dsQN <- normalizeQuantile(dsBG, path=normPath, name=name, typesToUpdate="pm", ..., verbose=verbose);

  verbose && exit(verbose);

  verbose && enter(verbose, "Computing expression summaries");
  
  if (summaryMethod=="rlm") {
    rmaPlm <- RmaPlm(dsQN);
    fitResult <- fit(rmaPlm, verbose=verbose);
  } else {
    throw("Only \"rlm\" is currently supported as a summary method");
  }

  verbose && exit(verbose);
  
  res <- this$fromFiles(path=getPath(rmaPlm), pattern="cel$", ...);

  res;
})


############################################################################
# HISTORY:
# 2006-10-10
# o Renamed rma and gcrma to rmaSummary and gcrmaSummary, to avoid clash
#   with existing functions.
# o Added gcrma() wrapper function.
# o Added rma() wrapper function.
############################################################################
