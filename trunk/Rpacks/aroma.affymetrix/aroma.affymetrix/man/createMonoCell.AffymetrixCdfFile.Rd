%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.createMonoCell.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{createMonoCell.AffymetrixCdfFile}
\alias{createMonoCell.AffymetrixCdfFile}
\alias{AffymetrixCdfFile.createMonoCell}
\alias{createMonoCell.AffymetrixCdfFile}
\alias{createMonoCell,AffymetrixCdfFile-method}


\title{Creates a mono-cell version of the CDF}

\description{
 Creates a mono-cell version of the CDF.
}

\usage{\method{createMonoCell}{AffymetrixCdfFile}(this, chipType=getChipType(this), suffix="monocell", sep="-", path="cdf", nbrOfCellsPerField=1, ..., ram=1, verbose=TRUE)}

\arguments{
  \item{chipType}{The chip type of the new CDF.}
  \item{suffix}{A suffix added to the chip type of the new CDF.}
  \item{sep}{A string separating the chip type and the suffix string.}
  \item{path}{The path where to store the new CDF file.}
  \item{nbrOfCellsPerField}{Number of cells per group field the new CDF
     should have.}
  \item{...}{Additional arguments passed to \code{\link[affxparser]{writeCdf}}.}
  \item{ram}{A \code{\link[base]{double}} saying if more or less units should be converted
     per chunk.  A smaller value uses less memory.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the monocell CDF as an \code{\link{AffymetrixCdfFile}} object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCdfFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
