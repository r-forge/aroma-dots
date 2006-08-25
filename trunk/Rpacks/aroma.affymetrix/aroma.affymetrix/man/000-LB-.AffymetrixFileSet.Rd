%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{[.AffymetrixFileSet}
\alias{[.AffymetrixFileSet}
\alias{AffymetrixFileSet.[}
\alias{[.AffymetrixFileSet}
\alias{[,AffymetrixFileSet-method}


\title{Extract a subset of the file set}

\description{
  Extract a subset of the file set.
  This is just a wrapper for \code{\link[aroma.affymetrix:extract.AffymetrixFileSet]{*extract}()}.
}

\usage{\method{[}{AffymetrixFileSet}(this, i=NULL, ...)}

\arguments{
 \item{i}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying the files to be returned.}
 \item{...}{Arguments passed to \code{\link[aroma.affymetrix:extract.AffymetrixFileSet]{*extract}()}.}
}

\value{
  Returns what \code{\link[aroma.affymetrix:extract.AffymetrixFileSet]{*extract}()} returns.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:extract.AffymetrixFileSet]{*extract}()}.
  \code{\link[aroma.affymetrix:000-LB--LB-.AffymetrixFileSet]{*[[}()}.
  For more information see \code{\link{AffymetrixFileSet}}.
}
\keyword{internal}
\keyword{methods}
