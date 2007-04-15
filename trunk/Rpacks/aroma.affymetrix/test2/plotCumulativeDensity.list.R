#########################################################################/**
# @set "class=list"
# @RdocMethod plotCumulativeDensity
# @alias plotCumulativeDensity.data.frame
# @alias plotCumulativeDensity.matrix
# @alias plotCumulativeDensity.numeric
#
# @title "Plots cumulative density distributions for a set of vector"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{A single of @list of @numeric @vectors, a @numeric @matrix,
#     or a @numeric @data.frame.}
#  \item{xlim,ylim}{@character @vector of length 2. The x and y limits.}
#  \item{xlab,ylab}{@character string for labels on x and y axis.}
#  \item{col}{The color(s) of the curves.}
#  \item{lty}{The types of curves.}
#  \item{lwd}{The width of curves.}
#  \item{...}{Additional arguments passed to @see "stats::ecdf",
#    @see "graphics::plot", and @see "graphics::lines".}
#  \item{add}{If @TRUE, the curves are plotted in the current plot,
#   otherwise a new is created.}
# }
#
# @author
#*/#########################################################################
setMethodS3("plotCumulativeDensity", "list", function(X, xlim=NULL, ylim=c(0,1), xlab=NULL, ylab="cumulative density", col=1:length(X), lty=NULL, lwd=NULL, spar=NA, ..., add=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'X':
  nbrOfSamples <- length(X);

  # Argument 'xlab':
  if (is.null(xlab))
    xlab <- substitute(X);

  # Argument 'col':
  if (is.null(col)) {
    col <- seq(length=nbrOfSamples);
  } else {
    col <- rep(col, length.out=nbrOfSamples);
  }

  # Argument 'lty':
  if (!is.null(lty))
    lty <- rep(lty, length.out=nbrOfSamples);

  # Argument 'lwd':
  if (!is.null(lwd))
    lwd <- rep(lwd, length.out=nbrOfSamples);

  if (is.null(xlim))
    xlim <- range(sapply(X, FUN=range, na.rm=TRUE));
  if (is.null(ylim))
    ylim <- c(0,1);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot the densities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (add == FALSE) {
    suppressWarnings({
      plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...);
    })
  }

  xs <- seq(from=xlim[1], to=xlim[2], length=1000);
  fits <- vector("list", nbrOfSamples);
  for(kk in seq(along=fits)) {
    x <- X[[kk]];
    fcn <- ecdf(x);
    fit <- list(
      x = xs,
      y = fcn(xs),
      ecdf = fcn
    );

    if (!is.na(spar)) {
      suppressWarnings({
        fit <- smooth.spline(fit, spar=spar);
      })
    }

#    fit$x <- c(xlim[1], fit$x, xlim[2]);
#    fit$y <- c(0, fit$y, 1);

    suppressWarnings({
      lines(fit, col=col[kk], lty=lty[kk], lwd=lwd[kk], ...);
    })
    fits[[kk]] <- fit;
  }

  invisible(fits);
}) # plotCumulativeDensity()



setMethodS3("plotCumulativeDensity", "data.frame", function(X, xlab=NULL, ...) {
  # Argument 'xlab':
  if (is.null(xlab))
    xlab <- substitute(X);
  plotCumulativeDensity(as.list(X), xlab=xlab, ...);
})



setMethodS3("plotCumulativeDensity", "matrix", function(X, xlab=NULL, ...) {
  # Argument 'xlab':
  if (is.null(xlab))
    xlab <- substitute(X);
  plotCumulativeDensity(as.data.frame(X), xlab=xlab, ...);
})


setMethodS3("plotCumulativeDensity", "numeric", function(X, xlab=NULL, ...) {
  # Argument 'xlab':
  if (is.null(xlab))
    xlab <- substitute(X);
  plotCumulativeDensity(list(X), xlab=xlab, ...);
})


##############################################################################
# HISTORY:
# 2007-04-14
# o Created from plotDensity().
##############################################################################

