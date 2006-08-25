###########################################################################/**
# @set "class=AffymetrixDataFile"
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
#      @seemethod "getProbes", @seemethod "getProbeIntensities" and 
#      @see "aroma.light::plotDensity.numeric".}
#   \item{xlab,ylab}{The labels on the x and the y axes.}
#   \item{log}{If @TRUE, the density of the log (base 2) values are 
#      used, otherwise the non-logged values.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
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
setMethodS3("plotDensity", "AffymetrixDataFile", function(this, probes=1/20, types=NULL, ..., xlab=NULL, ylab="density (integrates to one)", log=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'probes':

  # Argument 'xlab':
  if (is.null(xlab)) {
    if (log) {
      xlab <- expression(log[2](y));
    } else {
      xlab <- expression(y);
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  suppressWarnings({
    probes <- getProbes(this, probes=probes, types=types, verbose=verbose, ...);
  })
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot density
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Plotting the density");
  verbose && cat(verbose, "Array: ", getPathname(this));
  suppressWarnings({
    verbose && enter(verbose, "Loading probe intensities");
    y <- getProbeIntensities(this, probes=probes, ...);
    verbose && exit(verbose);
    if (log) {
      verbose && cat(verbose, "Taking the logarithm (base 2)");
      y <- log(y, base=2);
    }
    verbose && cat(verbose, "Plotting");
    plotDensity(y, xlab=xlab, ylab=ylab, ...);
  })
  verbose && exit(verbose);
})



############################################################################
# HISTORY:
# 2006-07-27
# o Added argument 'verbose' to plotDensity().
# 2006-05-29
# o Added Rdoc comments.
# 2006-05-16
# o Created.
############################################################################
