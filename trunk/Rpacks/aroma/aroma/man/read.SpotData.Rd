%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{read.SpotData}
\alias{read.SpotData}
\alias{SpotData.read}
\alias{read.SpotData}
\alias{read,SpotData-method}


\title{Reads several Spot files into a SpotData object}

\description{
  Reads several Spot files into a SpotData object.
}

\synopsis{read.SpotData(this, filename=NULL, path=NULL, pattern=NULL, verbose=TRUE, ...)}
\usage{SpotData$read(filename=NULL, path=NULL, pattern=NULL, verbose=TRUE, ...)}

\arguments{
  \item{filename}{A \code{\link[base]{vector}} of filenames. Either \code{pattern} or \code{filename} must be specified.}
  \item{path}{A string (or an optional \code{\link[base]{vector}} of paths if \code{filename} is specified) to the files.}
  \item{pattern}{A pattern string for matching filenames. Either \code{pattern} or \code{filename} must be specified.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, information will printed out during
                 the reading of the file.}
}

\value{Returns a \code{\link{SpotData}} object.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{\dontrun{# For an example see help(SpotData).}}

\seealso{
  For pattern formats see \code{\link[base]{list.files}}().
  For more information see \code{\link{SpotData}}.
}
\keyword{internal}
\keyword{methods}
