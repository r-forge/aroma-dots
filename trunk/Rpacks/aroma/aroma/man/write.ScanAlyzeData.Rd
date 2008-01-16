%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ScanAlyzeData.R
% 
% on Tue Jan 15 18:36:51 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{write.ScanAlyzeData}
\alias{write.ScanAlyzeData}
\alias{ScanAlyzeData.write}
\alias{write.ScanAlyzeData}
\alias{write,ScanAlyzeData-method}


\title{Write a ScanAlyze Results Data file}

\description{
 Writes the ScanAlyzeData object to a file with the ScanAlyze file format.
}

\usage{\method{write}{ScanAlyzeData}(this, filename, path=NULL, slide=NULL, overwrite=FALSE, row.names=FALSE, ..., verbose=FALSE)}

\arguments{
  \item{filename}{The filename of the file to be read.}
  \item{path}{The path to the file.}
  \item{slide}{An \code{\link[base]{integer}} specifying which slides to be written to file. Currently, only one slide at the time can be written.}
  \item{...}{Arguments passed to \code{write.table}.}
}

\value{Returns nothing.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  sa <- ScanAlyzeData$read("group4.dat", path=system.file("data-ex", package="aroma"))

  # Writes the ScanAlyzeData to a temporary file. Note that this
  # file won't be exactly the same since a few lines, specifying
  # for instance the creator of the file, will be added. The data,
  # however, will be the same.
  filename <- paste(tempfile("ScanAlyzeData"), ".dat", sep="")
  write(sa, filename)

  sa2 <- ScanAlyzeData$read(filename)
  print(equals(sa, sa2))  # TRUE

  unlink(filename)
}

\references{
  The data file in the example above was contributed by
  Michael Stadler, Bioinformatics Group, Swiss Institute for
  Experimental Cancer Research.
}

\seealso{
  To read one or more ScanAlyze Results files
  see \code{\link[aroma:read.ScanAlyzeData]{*read}()}.
  For more information see \code{\link{ScanAlyzeData}}.
}
\keyword{internal}
\keyword{methods}
