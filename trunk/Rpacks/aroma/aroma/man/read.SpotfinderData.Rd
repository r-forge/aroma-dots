%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotfinderData.R
% 
% on Tue Jan 15 18:37:29 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{read.SpotfinderData}
\alias{read.SpotfinderData}
\alias{SpotfinderData.read}
\alias{read.SpotfinderData}
\alias{read,SpotfinderData-method}


\title{Reads several Spotfinder/TAV files into a SpotfinderData object}

\description{
 Reads several Spotfinder/TAV files into a SpotfinderData object.
 The acronym TAV stands for \emph{TIGR Array Viewer}.
}

\synopsis{read.SpotfinderData(this, filename=NULL, path=NULL, pattern=NULL, verbose=TRUE, ...)}
\usage{SpotfinderData$read(filename=NULL, path=NULL, pattern=NULL, verbose=TRUE, ...)}

\arguments{
  \item{filename}{A \code{\link[base]{vector}} of filenames. Either \code{pattern} or \code{filename} must be specified.}
  \item{path}{A string (or an optional \code{\link[base]{vector}} of paths if \code{filename} is specified) to the files.}
  \item{pattern}{A pattern string for matching filenames. Either \code{pattern} or \code{filename} must be specified.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, information will printed out during
                 the reading of the file.}
}

\value{Returns a \code{\link{SpotfinderData}} object.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{\dontrun{# For an example see help(SpotfinderData).}}

\seealso{
  For pattern formats see \code{\link[base]{list.files}}().
  For more information see \code{\link{SpotfinderData}}.
}
\keyword{internal}
\keyword{methods}
