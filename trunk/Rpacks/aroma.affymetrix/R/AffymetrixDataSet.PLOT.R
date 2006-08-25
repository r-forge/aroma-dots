###########################################################################/**
# @set "class=AffymetrixDataSet"
# @RdocMethod plotDensity
#
# @title "Plots the densities of all samples"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{probes}{The subset of probes to include.}
#   \item{types}{The type of probes to include.}
#   \item{...}{Additional arguments passed to 
#      @see "plotDensity.AffymetrixDataFile".}
#   \item{col}{A @vector of colors for each of the arrays.}
#   \item{lty}{A @vector of line types for each of the arrays.}
#   \item{lwd}{A @vector of line widths for each of the arrays.}
#   \item{add}{If @FALSE, a new plot is created, otherwise the generated
#     graphics is added to the current plot.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("plotDensity", "AffymetrixDataSet", function(this, probes=1/20, types=NULL, ..., col=seq(this), lty=NULL, lwd=NULL, add=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'probes':

  nbrOfArrays <- nbrOfArrays(this);
  
  # Argument 'col':
  if (is.null(col)) {
    col <- seq(length=nbrOfArrays);
  } else {
    col <- rep(col, length.out=nbrOfArrays);
  }

  # Argument 'lty':
  if (!is.null(lty))
    lty <- rep(lty, length.out=nbrOfArrays);

  # Argument 'lwd':
  if (!is.null(lwd))
    lwd <- rep(lwd, length.out=nbrOfArrays);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  suppressWarnings({
    probes <- getProbes(this, probes=probes, types=types, ...);
  })


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot densities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  dataFiles <- this$dataFiles;
  for (cc in seq(length=nbrOfArrays)) {
    add <- add || (cc > 1);
    plotDensity(dataFiles[[cc]], probes=probes, ..., col=col[cc], lty=lty[cc], lwd=lwd[cc], add=add);
  }
})



############################################################################
# HISTORY:
# 2006-05-16
# o Created.
############################################################################
