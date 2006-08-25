%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ParameterSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{[[.ParameterSet}
\alias{[[.ParameterSet}
\alias{ParameterSet.[[}
\alias{[[.ParameterSet}
\alias{[[,ParameterSet-method}


\title{Gets a subset of the files in the parameter set}

\description{
  Gets a subset of the files in the parameter set.
}

\usage{\method{[[}{ParameterSet}(this, i=NULL, drop=TRUE)}

\arguments{
 \item{i}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying files to be returned.}
 \item{drop}{If \code{\link[base:logical]{TRUE}} and only one array is returned, the data file is
   return directly, otherwise as an element in a list.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{list}} of \code{\link{ParameterFile}} objects, or a single
  \code{\link{ParameterFile}} object if \code{drop} is \code{\link[base:logical]{TRUE}} and only
  one data file was selected.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:000-LB-.ParameterSet]{*[}()}.
  For more information see \code{\link{ParameterSet}}.
}
\keyword{internal}
\keyword{methods}
