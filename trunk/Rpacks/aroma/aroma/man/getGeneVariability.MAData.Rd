%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getGeneVariability.MAData}
\alias{getGeneVariability.MAData}
\alias{MAData.getGeneVariability}
\alias{getGeneVariability.MAData}
\alias{getGeneVariability,MAData-method}


\title{Gets the genewise variability}

\usage{\method{getGeneVariability}{MAData}(this, robust=TRUE, force=FALSE, slides=NULL, ...)}

\arguments{
  \item{robust}{If \code{\link[base:logical]{TRUE}} the median absolute deviation (MAD) of
    the residuals will be calculated, otherwise the sample standard
    deviation will be calculated.}
  \item{force}{If \code{\link[base:logical]{FALSE}} and if cached gene variability values
    exists they will be used, otherwise the gene variability will be
    (re-)calculated.}
  \item{slides}{The slides which should be included in the calculations.
    If \code{\link[base]{NULL}}, all slides are included.}
}

\description{
  Calculates the genewise variability. What spots belongs to which genes
  is defined by the layout of the slides, where all slides are assumed to
  have the same layout. See class \code{\link{Layout}} for more
  information about genes. For speed improvement, the gene variability
  will be cached for future queries. To override the cache, use
  \code{force=TRUE}.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw, bgSubtract=TRUE)
  ma.norm <- clone(ma)
  normalizeWithinSlide(ma.norm, method="s")
  normalizeAcrossSlides(ma.norm)
  var <- getGeneVariability(ma)
  var.norm <- getGeneVariability(ma.norm)
  var.diff <- var - var.norm

  # Statistics
  print(summary(var))
  print(summary(var.norm))
  print(summary(var.diff))
  cat("Number of 'improved' genes :", sum(var.diff > 0, na.rm=TRUE), "\n")
  cat("Number of 'worsened' genes :", sum(var.diff < 0, na.rm=TRUE), "\n")
  cat("Number of 'unchanged' genes:", sum(var.diff == 0, na.rm=TRUE), "\n")

  # Plots
  xlim <- quantile(c(var, var.norm), probs=c(0,0.999))
  subplots(4)
  hist(var, nclass=100, xlim=xlim, cex.main=0.7,
            main="Gene variability before normalization");
  hist(var.norm, nclass=100, xlim=xlim, cex.main=0.7,
                 main="Gene variability after normalization");
  hist(var.diff, nclass=50, cex.main=0.7, main="Improvements");
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  To calculate the mean (or any other quantile) of the genewise
  variabilities see \code{\link[aroma:getMOR.MAData]{*getMOR}()}.
  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
