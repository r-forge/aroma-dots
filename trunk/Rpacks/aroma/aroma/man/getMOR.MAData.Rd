%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getMOR.MAData}
\alias{getMOR.MAData}
\alias{MAData.getMOR}
\alias{getMOR.MAData}
\alias{getMOR,MAData-method}


\title{Gets the Measure of Reproducibility}

\usage{\method{getMOR}{MAData}(this, robust=TRUE, probs=NULL, force=FALSE, slides=NULL, ...)}

\arguments{
  \item{robust}{If \code{\link[base:logical]{TRUE}} the median absolute deviation (MAD) will
    be used for calculating the genewise residuals, otherwise the sample
    standard deviation will be used.}
  \item{probs}{The quantiles of the genewise variabilities returned. If
    \code{-0.5} (or \code{\link[base]{NULL}}) the mean is returned.}
  \item{force}{If \code{\link[base:logical]{FALSE}} and if cached gene variability values exists
    they will be used, otherwise the gene variability will be
    (re-)calculated.}
  \item{slides}{The slides which should be included in the calculations.
    If \code{\link[base]{NULL}}, all slides are included.}
}

\description{
  Calculates the Measure of Reproducibility (MOR) [1], which is by default
  the (scalar) mean value of all genewise robust variability measures
  given by \code{getGeneVariability()}, either the standard deviation or
  the median absolute deviation (MAD). Any quantile(s) of these can be
  returned by setting argument \code{probs}.
}

\value{
  Returns a \code{\link[base]{vector}} of the quantiles of the genewise variabilities asked
  for. By default the mean value of all genewise variabilities (MOR) is
  returned.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw, bgSubtract=TRUE)
  ma.norm <- clone(ma)
  normalizeWithinSlide(ma.norm, method="s")
  normalizeAcrossSlides(ma.norm)
  mor <- getMOR(ma)
  mor.norm <- getMOR(ma.norm)

  cat("MOR for non-normalized data:", mor, "\n")
  cat("MOR for normalized data:", mor.norm, "\n")
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
 [1] Henrik Bengtsson, Plate Effects in cDNA microarray data,
     Matemathical Statistics, Centre for Matematical Sciences,
     Lund University, Sweden. Manuscript, 2002.
}

\seealso{
  \code{\link[aroma:getVariability.MAData]{*getVariability}()}
  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
