%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{as.data.frame.MicroarrayData}
\alias{as.data.frame.MicroarrayData}
\alias{MicroarrayData.as.data.frame}
\alias{as.data.frame.MicroarrayData}
\alias{as.data.frame,MicroarrayData-method}


\title{Converts the object to a data frame}

\description{
  Converts the object to a data frame.
}

\usage{\method{as.data.frame}{MicroarrayData}(x, ..., optional=TRUE)}

\arguments{
 \item{...}{Not used.}
 \item{optional}{If \code{\link[base:logical]{TRUE}}, optional fields are also returned.}
}

\value{
  Return a \code{\link[base]{data.frame}}.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  df <- as.data.frame(raw)
}

\seealso{
 \code{\link[aroma:extract.MicroarrayData]{*extract}()},
 \code{\link[base]{data.frame}}().
  For more information see \code{\link{MicroarrayData}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
