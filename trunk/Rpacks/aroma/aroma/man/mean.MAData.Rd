%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.R
% 
% on Tue Jan 15 18:33:22 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{mean.MAData}
\alias{mean.MAData}
\alias{MAData.mean}
\alias{mean.MAData}
\alias{mean,MAData-method}


\title{Average Mean for microarray data}

\description{
  Computes the average mean across all slides for each spot.
}

\usage{\method{mean}{MAData}(x, inf.rm=TRUE, na.rm=TRUE, ...)}

\arguments{
 \item{inf.rm}{a logical value indicating whether \code{\link[base:is.finite]{Inf}} values should be
        stripped before the computation proceeds.}
 \item{na.rm}{a logical value indicating whether \code{\link[base]{NA}} values should be
        stripped before the computation proceeds.}
 \item{...}{other arguments to \code{mean}, e.g. \code{trim}.}
}

\value{Returns a \code{\link[base]{matrix}} with the columns \code{meanM} and \code{meanA}.}

\examples{
  # The option 'dataset' is used to annotate plots.
  options(dataset="sma:MouseArray")

  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw, bgSubtract=TRUE)
  mean <- mean(ma)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:var.MAData]{*var}()}.
  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
