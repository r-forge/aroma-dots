%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ScanAlyzeData.R
% 
% on Tue Jan 15 18:36:51 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{read.ScanAlyzeData}
\alias{read.ScanAlyzeData}
\alias{ScanAlyzeData.read}
\alias{read.ScanAlyzeData}
\alias{read,ScanAlyzeData-method}


\title{Reads one or several ScanAlyze files into a ScanAlyzeData object}

\description{
 Reads one or several ScanAlyze files into a ScanAlyzeData object.
}

\synopsis{read.ScanAlyzeData(static, filename=NULL, path=NULL, pattern=NULL, verbose=FALSE, ...)}
\usage{ScanAlyzeData$read(filename=NULL, path=NULL, pattern=NULL, verbose=FALSE, ...)}

\arguments{
  \item{filename}{A \code{\link[base]{vector}} of filenames. Either \code{pattern} or \code{filename} must be specified.}
  \item{path}{A string (or an optional \code{\link[base]{vector}} of paths if \code{filename} is specified) to the files.}
  \item{pattern}{A pattern string for matching filenames. Either \code{pattern} or \code{filename} must be specified.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, information will printed out during
                 the reading of the file.}
}

\value{Returns a \code{\link{ScanAlyzeData}} object.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  sa <- ScanAlyzeData$read(pattern="group.*.dat", path=system.file("data-ex", package="aroma"))
  raw <- getRawData(sa)
}

\references{
  The two data files in the example above was contributed by
  Michael Stadler, Bioinformatics Group, Swiss Institute for
  Experimental Cancer Research.
}

\seealso{
  For pattern formats see \code{\link[base]{list.files}}().
  For more information see \code{\link{ScanAlyzeData}}.
}
\keyword{internal}
\keyword{methods}
