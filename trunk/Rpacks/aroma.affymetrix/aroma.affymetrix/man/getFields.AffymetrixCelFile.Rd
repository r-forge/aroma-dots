%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getFields.AffymetrixCelFile}
\alias{getFields.AffymetrixCelFile}
\alias{AffymetrixCelFile.getFields}
\alias{getFields.AffymetrixCelFile}
\alias{getFields,AffymetrixCelFile-method}


\title{Gets all or a subset of the fields in a CEL file}

\description{
 Gets all or a subset of the fields in a CEL file for all or a subset of the cells.
}

\usage{\method{getFields}{AffymetrixCelFile}(this, indices=NULL, fields=c("xy", "intensities", "stdvs", "pixels"), ..., verbose=FALSE)}

\arguments{
  \item{indices}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of cell indices.  If \code{\link[base]{NULL}}, all cells
    are considered.}
  \item{fields}{Names of fields to be retrieved.}
  \item{...}{Additional arguments passed to \code{\link[affxparser]{readCel}}.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{data.frame}} of the fields requested.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
