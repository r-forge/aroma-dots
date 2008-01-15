%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readHeader.GenePixData}
\alias{readHeader.GenePixData}
\alias{GenePixData.readHeader}
\alias{readHeader.GenePixData}
\alias{readHeader,GenePixData-method}


\title{Reads a GenePix Results (GPR) file header}

\description{
  Reads a GenePix Results (GPR) file header. This method will parse the header of a GPR file and
  extract all header information including the column names of the
  preceeding data section.
}

\synopsis{readHeader.GenePixData(static, filename, path=NULL, verbose=FALSE, ...)}
\usage{GenePixData$readHeader(filename, path=NULL, verbose=FALSE, ...)}

\arguments{
  \item{filename}{The filename of the file to be read.}
  \item{path}{The path to the file.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, information will printed out during
                 the reading of the file.}
}

\value{Returns a \code{\link[base]{list}} of the header fields, but also the names and
  the "guessed" classes of the columns to be read.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\details{
  This method is used internally by the \code{\link[aroma:read.GenePixData]{*read}()} method.

  There are many difference version of the GenePix image analysis
  software and almost equally number of different versions of the
  GenePix Results file format. All example data available at
  \url{http://www.axon.com/gn_GPR_Format_History.html} have been
  tested and can successfully be read using this method.
}

\references{
  GenePix File Formats:
  \url{http://www.axon.com/GN_GenePix_File_Formats.html},
  \url{http://www.axon.com/gn_GPR_Format_History.html}
}

\examples{
  header <- GenePixData$readHeader("gpr123.gpr", path=system.file("data-ex", package="aroma"))
  str(header)
}

\seealso{
  To read a GenePix Results file including the data section
  see \code{\link[aroma:read.GenePixData]{*read}()}.
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
