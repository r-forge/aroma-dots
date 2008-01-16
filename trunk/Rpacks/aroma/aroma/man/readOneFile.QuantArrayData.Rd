%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QuantArrayData.R
% 
% on Tue Jan 15 18:35:00 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readOneFile.QuantArrayData}
\alias{readOneFile.QuantArrayData}
\alias{QuantArrayData.readOneFile}
\alias{readOneFile.QuantArrayData}
\alias{readOneFile,QuantArrayData-method}


\title{Reads a QuantArray result file}

\description{
 Reads a QuantArray result file.
 Currently QuantArray v2 and v3 are supported. If the version is not
 recognized it will try to make a best guess and read it in anyway.
 However, it can not be guaranteed that all fields will be of the correct
 data type.

 This method also reads Unicoded QuantArray files. If the file is
 Unicoded it is first translated into a temporary ASCII file, which is
 then read. The translation from Unicode to ASCII is done by "brute force",
 i.e. by excluding the highbyte and only keeping the lowbytes. This means
 that some characters will be incorrectly translated. If that happens, a
 warning will be given, otherwise not.
}

\synopsis{readOneFile.QuantArrayData(this, filename, path=NULL, verbose=FALSE, ...)}
\usage{QuantArrayData$readOneFile(filename, path=NULL, verbose=FALSE, ...)}

\arguments{
  \item{filename}{The filename of the file to be read.}
  \item{path}{The path to the file.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, information will printed out during
                 the reading of the file.}
}

\value{Returns a \code{\link{QuantArrayData}} object.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
 \dontrun{
  # At the moment there is no QuantArray sample files in the package...
  # qa <- QuantArrayData$read("quantarray123.txt", path=system.file("data-ex", package="aroma"))
  # ...will use a GenePix sample file instead.
  qa <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))

  raw <- getRawData(qa)
 }
}

\seealso{
  For more information see \code{\link{QuantArrayData}}.
}
\keyword{internal}
\keyword{methods}
