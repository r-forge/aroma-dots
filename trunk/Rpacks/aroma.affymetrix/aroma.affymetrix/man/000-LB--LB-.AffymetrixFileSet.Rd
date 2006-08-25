%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{[[.AffymetrixFileSet}
\alias{[[.AffymetrixFileSet}
\alias{AffymetrixFileSet.[[}
\alias{[[.AffymetrixFileSet}
\alias{[[,AffymetrixFileSet-method}


\title{Gets a subset of the files in the file set}

\description{
  Gets a subset of the files in the file set.
}

\usage{\method{[[}{AffymetrixFileSet}(this, i=NULL, drop=TRUE)}

\arguments{
 \item{i}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying which files to be returned.}
 \item{drop}{If \code{\link[base:logical]{TRUE}} and only one file is returned, the file is
   return directly, otherwise as an element in a list.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{list}} of \code{\link{AffymetrixFile}} objects, or a single
  \code{\link{AffymetrixFile}} object if \code{drop} is \code{\link[base:logical]{TRUE}} and only
  one file was selected.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:000-LB-.AffymetrixFileSet]{*[}()}.
  For more information see \code{\link{AffymetrixFileSet}}.
}
\keyword{internal}
\keyword{methods}
