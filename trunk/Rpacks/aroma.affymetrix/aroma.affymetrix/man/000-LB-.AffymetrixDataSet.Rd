%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{[.AffymetrixDataSet}
\alias{[.AffymetrixDataSet}
\alias{AffymetrixDataSet.[}
\alias{[.AffymetrixDataSet}
\alias{[,AffymetrixDataSet-method}


\title{Extract a subset of the data set}

\description{
  Extract a subset of the data set.
  This is just a wrapper for \code{\link[aroma.affymetrix:extract.AffymetrixDataSet]{*extract}()}.
}

\usage{\method{[}{AffymetrixDataSet}(this, i=NULL, ...)}

\arguments{
 \item{i}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying the arrays to be returned.}
 \item{...}{Arguments passed to \code{\link[aroma.affymetrix:extract.AffymetrixDataSet]{*extract}()}.}
}

\value{
  Returns what \code{\link[aroma.affymetrix:extract.AffymetrixDataSet]{*extract}()} returns.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:extract.AffymetrixDataSet]{*extract}()}.
  \code{\link[aroma.affymetrix:000-LB--LB-.AffymetrixDataSet]{*[[}()}.
  For more information see \code{\link{AffymetrixDataSet}}.
}
\keyword{internal}
\keyword{methods}
