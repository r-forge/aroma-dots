%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{write.GenePixData}
\alias{write.GenePixData}
\alias{GenePixData.write}
\alias{write.GenePixData}
\alias{write,GenePixData-method}


\title{Write a GenePix Results Data file}

\description{
 Writes the GenePixData object to a file with the GenePix file format
 GPR (GenePix Results Format).
}

\usage{\method{write}{GenePixData}(this, filename, path=NULL, slide=1, overwrite=FALSE, ..., verbose=FALSE)}

\arguments{
  \item{filename}{The filename of the file to be written.}
  \item{path}{The path to the file.}
  \item{slide}{An integer specifying which slide to be written to file.
    Currently, only one slide at the time can be written.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already an existing file will be overwritten.
      Otherwise an Exception will be thrown.}
  \item{...}{Arguments passed to \code{write.table}.}
}

\value{Returns nothing.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  GenePix File Formats,
  \url{http://www.axon.com/GN_GenePix_File_Formats.html}
}

\examples{
  gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))

  # Writes the GenePix Results Data to a file named "temp.gpr". Note
  # that this file won't be exactly the same since a few lines,
  # specifying for instance the creator of the file, will be added.
  write(gpr, "temp.gpr", overwrite=TRUE)
  file.show("temp.gpr")
}

\seealso{
  To read one or more GenePix Results files
  see \code{\link[aroma:read.GenePixData]{*read}()}.
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
