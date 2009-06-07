%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{var.MAData}
\alias{var.MAData}
\alias{MAData.var}
\alias{var.MAData}
\alias{var,MAData-method}


\title{Variance for microarray data}

\description{
  Computes the variance across all slides for each spot.
}

\usage{\method{var}{MAData}(this, inf.rm=TRUE, na.rm=TRUE, ...)}

\arguments{
 \item{inf.rm}{a logical value indicating whether \code{\link[base:is.finite]{Inf}} values should be
        stripped before the computation proceeds.}
 \item{na.rm}{a logical value indicating whether \code{\link[base]{NA}} values should be
        stripped before the computation proceeds.}
 \item{...}{other arguments to \code{var}, e.g. \code{trim}.}
}

\value{Returns a \code{\link[base]{matrix}} with the columns \code{varM} and \code{varA}.}

\examples{
  # The option 'dataset' is used to annotate plots.
  options(dataset="sma:MouseArray")

  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw, bgSubtract=TRUE)
  var <- ma$var()
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:mean.MAData]{*mean}()}.
  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
