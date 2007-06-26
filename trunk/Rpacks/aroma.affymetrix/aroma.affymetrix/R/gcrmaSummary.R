###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod gcrmaSummary
#
# @title "Calculates the GCRMA expression summary"
#
# \description{
#  @get "title".
#
# Applies GC-based background correction, followed by quantile normalisation
# and then a summary by either median polish or robust regression.
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
#   \item{type}{The type of background correction.  Currently accepted types
#       are "fullmodel" (the default, uses MMs) and "affinities" (uses
#       probe sequence only).}
#   \item{indicesNegativeControl}{Locations of any negative control
#       probes (e.g., the anti-genomic controls on the human exon array).
#       If @NULL and type=="affinities", MMs are used as the negative
#       controls.}
#   \item{opticalAdjust}{If @TRUE, apply correction for optical effect,
#       as in @see "gcrma::bg.adjust.optical".}
#   \item{gsbAdjust}{Should we adjust for specific binding (defaults to
#        @TRUE)?}
#   \item{k}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
#   \item{rho}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
#   \item{stretch}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
#   \item{fast}{If @TRUE, an ad hoc transformation of the PM is performed
#       (\code{gcrma::gcrma.bg.transformation.fast}).}
#   \item{overwrite}{If @TRUE, already adjusted arrays are overwritten,
#     unless skipped, otherwise an error is thrown.}
#   \item{skip}{If @TRUE, the array is not normalized if it already exists.}
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
#  @see "gcrma::bg.adjust.gcrma"
#  @seeclass
# }
#*/###########################################################################
setMethodS3("gcrmaSummary", "AffymetrixCelSet", function(this, path=NULL, name="gcrma", bgPath=NULL, normPath=NULL, summaryMethod="rlm", probePath=NULL, affinities=NULL, type="fullmodel",  indicesNegativeControl=NULL, opticalAdjust=TRUE, gsbAdjust=TRUE, k=6 * fast + 0.5 * (1 - fast), rho=0.7, stretch=1.15*fast + (1-fast), fast=TRUE, ..., verbose=FALSE) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  if (is.null(path)) {
    # Path structure: /gcrma/<data set name>/chip_data/<chip type>/
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

  dsBG <- bgAdjustGcrma(this, path=bgPath, name=name, probePath=probePath, affinities=affinities, type=type,  indicesNegativeControl=indicesNegativeControl, opticalAdjust=opticalAdjust, gsbAdjust=gsbAdjust, k=6 * fast + 0.5 * (1 - fast), rho=0.7, stretch=1.15*fast + (1-fast), fast=fast, ..., verbose=verbose);

  verbose && exit(verbose);

  # normalisation step
  
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
