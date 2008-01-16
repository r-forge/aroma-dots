%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotPosition.R
% 
% on Tue Jan 15 18:37:36 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{write.SpotPosition}
\alias{write.SpotPosition}
\alias{SpotPosition.write}
\alias{write.SpotPosition}
\alias{write,SpotPosition-method}


\title{Write a SpotPosition object to file}

\description{
  Write a SpotPosition object to file. By default, if not overridden
  by a method in a subclass, it writes the data to a tab-delimited file.
  Note that subclasses like GenePixData, ScanAlyzeData and SpotData do
  write files in their special file formats. To force such object of such
  classes to be written as tab-delimited file, do
  \code{write.MicroarrayData(object, ...)} instead.
}

\usage{\method{write}{SpotPosition}(this, filename, path=NULL, slides=NULL, overwrite=FALSE, row.names=FALSE, sep="	", ..., verbose=FALSE)}

\arguments{
  \item{filename}{The filename of the GPR file to be written.}
  \item{path}{The path to the GPR file.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, an existing file is overwritten.
    Otherwise an exception is thrown.}
  \item{row.names}{If \code{\link[base:logical]{TRUE}}, row names are written, otherwise not.}
  \item{sep}{The separator between the cells.}
  \item{...}{Other arguments accepted by subclasses or which are passed
    to \code{write.table}.}
}

\value{Returns nothing.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:read.SpotPosition]{*read}()}.
  For more information see \code{\link{SpotPosition}}.
}
\keyword{internal}
\keyword{methods}
