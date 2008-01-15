%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{read.GenePixData}
\alias{read.GenePixData}
\alias{GenePixData.read}
\alias{read.GenePixData}
\alias{read,GenePixData-method}


\title{Reads one or several GenePix files into one GenePixData object}

\description{
 Reads one or several GenePix files into one GenePixData object.
}

\synopsis{read.GenePixData(this, filename=NULL, path=NULL, pattern=NULL, verbose=FALSE, ...)}
\usage{GenePixData$read(filename=NULL, path=NULL, pattern=NULL, verbose=FALSE, ...)}

\arguments{
  \item{filename}{A vector of filenames. Either \code{pattern} or \code{filename} must be specified.}
  \item{path}{A string (or an optional vector of paths if \code{filename} is specified) to the files.}
  \item{pattern}{A pattern string for matching filenames. Either \code{pattern} or \code{filename} must be specified.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, information will printed out during
                 the reading of the file.}
}

\value{Returns a \code{\link{GenePixData}} object.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  # Reads one GenePix data files
  gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))

  # Reads two GenePix data files
  filenames <- c("gpr123.gpr", "gpr123.gpr");
  gpr <- GenePixData$read(filenames, path=system.file("data-ex", package="aroma"))

  # Reads all GenePix data files
  pattern <- c(".*\\.gpr");
  gpr <- GenePixData$read(pattern=pattern, path=system.file("data-ex", package="aroma"))
}

\seealso{
  To write a slide to a GenePix Results file see \code{\link[aroma:write.GenePixData]{*write}()}.
  For pattern formats see \code{\link[base]{list.files}}().
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
