###########################################################################/**
# @set "class=AffymetrixCelSet"
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
#   \item{subset}{The subset of probes to include.}
#   \item{types}{The type of probes to include.}
#   \item{...}{Additional arguments passed to 
#      @see "plotDensity.AffymetrixCelFile".}
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
setMethodS3("plotDensity", "AffymetrixCelSet", function(this, subset=1/2, types=NULL, ..., col=seq(this), lty=NULL, lwd=NULL, annotate=TRUE, add=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'subset':

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

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);
  verbose && enter(verbose, "Identifying subset of probes");
  suppressWarnings({
    subset <- identifyCells(cdf, indices=subset, types=types, ...,
                                                    verbose=less(verbose));
  })
  verbose && exit(verbose);
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot densities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  for (cc in seq(length=nbrOfArrays)) {
    df <- getFile(this, cc);
    add <- add || (cc > 1);
    plotDensity(df, subset=subset, ..., col=col[cc], lty=lty[cc], 
              lwd=lwd[cc], annotate=FALSE, add=add, verbose=less(verbose));
    if (annotate) {
      stextChipType(getCdf(this));
      stextSize(df, size=length(subset));
      annotate <- FALSE;
    }
  }
})



############################################################################
# HISTORY:
# 2006-05-16
# o Created.
############################################################################
