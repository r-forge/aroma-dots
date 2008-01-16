%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QuantArrayData.R
% 
% on Tue Jan 15 18:35:00 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{read.QuantArrayData}
\alias{read.QuantArrayData}
\alias{QuantArrayData.read}
\alias{read.QuantArrayData}
\alias{read,QuantArrayData-method}


\title{Reads one or several QuantArray files into one QuantArrayData object}

\description{
 Reads one or several QuantArray files into one QuantArrayData object.
 Currently QuantArray v2 and v3 are supported. If the version is not
 recognized it will try to make a best guess and read it in anyway.
 However, it can not be guaranteed that all fields will be of the correct
 data type.

 This method also reads Unicoded QuantArray files. If the file is
 Unicoded it is first translated into a temporary ASCII file, which is
 then read. The translation from Unicode to ASCII is done by "brute force",
 i.e. by excluding the hi-bytes and only keeping the lo-bytes. This means
 that some characters will be incorrectly translated. If that happens, a
 warning will be given, otherwise not.
}

\synopsis{read.QuantArrayData(this, filename=NULL, path=NULL, pattern=NULL, verbose=FALSE, ...)}
\usage{QuantArrayData$read(filename=NULL, path=NULL, pattern=NULL, verbose=FALSE, ...)}

\arguments{
  \item{filename}{A \code{\link[base]{vector}} of filenames. Either \code{pattern} or \code{filename} must be specified.}
  \item{path}{A string (or an optional \code{\link[base]{vector}} of paths if \code{filename} is specified) to the files.}
  \item{pattern}{A pattern string for matching filenames. Either \code{pattern} or \code{filename} must be specified.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, information will printed out during
                 the reading of the file.}
}

\value{Returns a \code{\link{QuantArrayData}} object.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
 \dontrun{
  # Read two QuantArray data files
  # filenames <- c("quantarray123.gpr", "quantarray123.gpr");
  # qa <- QuantArrayData$read(filenames, path=system.file("data-ex", package="aroma"))

  # Do not have any QuantArray sample files, will examplify using
  # the GenePix sample files instead...
  filenames <- c("gpr123.gpr", "gpr123.gpr");
  qa <- GenePixData$read(filenames, path=system.file("data-ex", package="aroma"))

  # Create a RawData object from this QuantArrayData objects.
  raw <- getRawData(qa)
 }
}

\seealso{
  To write a slide to a QuantArray Results file
  see \code{\link[aroma:write.QuantArrayData]{*write}()}.
  For pattern formats see \code{\link[base]{list.files}}().
  For more information see \code{\link{QuantArrayData}}.
}
\keyword{internal}
\keyword{methods}
