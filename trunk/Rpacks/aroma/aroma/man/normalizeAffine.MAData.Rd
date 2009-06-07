%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeAffine.MAData}
\alias{normalizeAffine.MAData}
\alias{MAData.normalizeAffine}
\alias{normalizeAffine.MAData}
\alias{normalizeAffine,MAData-method}


\title{Affine normalization based on non-logged data}

\description{
  Affine normalization based on non-logged data.
  For details, see \code{\link{normalizeAffine.RGData}}().
}

\usage{\method{normalizeAffine}{MAData}(this, ...)}

\arguments{
 \item{...}{Arguments passed to \code{\link{normalizeAffine.RGData}}().}
}

\examples{
  # The option 'dataset' is used to annotate plots.
  options(dataset="sma:MouseArray")

  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)

  # Clone the data to get one non-normalized and one normalized data set.
  maCurveFit <- clone(ma)
  maAffine <- clone(ma)

  # Normalize the data within slides using scaled print-tip normalization.
  normalizeCurveFit(maCurveFit, groupBy="printtip", method="lowess")
  normalizeAffine(maAffine, groupBy="printtip")

  # Plot data before and after normalization.
  subplots(9, nrow=3)
  for (obj in list(ma, maCurveFit, maAffine)) {
    # Plot M vs A and M spatially for array 1.
    plot(obj)
    plotSpatial(obj)
    # Plot the densities of A for *all* arrays.
    plotDensity(obj, what="A", xlim=c(4,16))
  }
}

\seealso{
  For more information see \code{\link{MAData}}.
}

\author{
  Henrik Bengtsson (\url{http://www.braju.com/R/}).
}
\keyword{internal}
\keyword{methods}
