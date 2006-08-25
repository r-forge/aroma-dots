%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{[[.AffymetrixDataSet}
\alias{[[.AffymetrixDataSet}
\alias{AffymetrixDataSet.[[}
\alias{[[.AffymetrixDataSet}
\alias{[[,AffymetrixDataSet-method}


\title{Gets a subset of the data files in the data set}

\description{
  Gets a subset of the data files in the data set.
}

\usage{\method{[[}{AffymetrixDataSet}(this, i=NULL, drop=TRUE)}

\arguments{
 \item{i}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying data files to be returned.}
 \item{drop}{If \code{\link[base:logical]{TRUE}} and only one array is returned, the data file is
   return directly, otherwise as an element in a list.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{list}} of \code{\link{AffymetrixDataFile}} objects, or a single
  \code{\link{AffymetrixDataFile}} object if \code{drop} is \code{\link[base:logical]{TRUE}} and only
  one data file was selected.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:000-LB-.AffymetrixDataSet]{*[}()}.
  For more information see \code{\link{AffymetrixDataSet}}.
}
\keyword{internal}
\keyword{methods}
