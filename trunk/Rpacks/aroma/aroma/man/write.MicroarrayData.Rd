%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.IO.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{write.MicroarrayData}
\alias{write.MicroarrayData}
\alias{MicroarrayData.write}
\alias{write.MicroarrayData}
\alias{write,MicroarrayData-method}


\title{Write a MicroarrayData object to file}

\description{
  Write a MicroarrayData object to file. By default, if not overridden
  by a method in a subclass, it writes the data to a tab-delimited file.
  Note that subclasses like GenePixData, ScanAlyzeData and SpotData do
  write files in their special file formats. To force such object of such
  classes to be written as tab-delimited file, do
  \code{write.MicroarrayData(object, ...)} instead.
}

\usage{\method{write}{MicroarrayData}(this, filename, path=NULL, slides=NULL, overwrite=FALSE, row.names=FALSE, sep="	", ..., verbose=FALSE)}

\arguments{
  \item{filename}{The filename of the GPR file to be written.}
  \item{path}{The path to the GPR file.}
  \item{slides}{The slides to be written. If \code{\link[base]{NULL}}, all slides are considered.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, an existing file is overwritten.
    Otherwise an exception is thrown.}
  \item{row.names}{If \code{\link[base:logical]{TRUE}}, row names are written, otherwise not.}
  \item{sep}{The separator between the cells.}
  \item{...}{Other arguments accepted by subclasses or which are passed
    to \code{write.table}.}
}

\value{Returns nothing.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  # Loads the file 'gpr123.gpr' located in the data directory:
  gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))

  # Writes the GenePix Results Data to a file named "temp.gpr". Note
  # that this file won't be exactly the same since a few lines,
  # specifying for instance the creator of the file, will be added.
  write(gpr, "temp.gpr", overwrite=TRUE)

  # Extracts the raw data from the gpr object and saves it to file.
  raw <- getRawData(gpr)
  write(raw, "temp.raw", overwrite=TRUE)
  file.show("temp.raw")
}

\seealso{
  To read one or more MicroarrayData files at once see \code{\link[aroma:read.MicroarrayData]{*read}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
