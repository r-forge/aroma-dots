%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ParameterSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{[.ParameterSet}
\alias{[.ParameterSet}
\alias{ParameterSet.[}
\alias{[.ParameterSet}
\alias{[,ParameterSet-method}


\title{Extract a subset of the parameter set}

\description{
  Extract a subset of the parameter set.
  This is just a wrapper for \code{\link[aroma.affymetrix:extract.ParameterSet]{*extract}()}.
}

\usage{\method{[}{ParameterSet}(this, i=NULL, ...)}

\arguments{
 \item{i}{An \code{\link[base]{integer}} index \code{\link[base]{vector}} specifying the arrays to be returned.}
 \item{...}{Arguments passed to \code{\link[aroma.affymetrix:extract.ParameterSet]{*extract}()}.}
}

\value{
  Returns what \code{\link[aroma.affymetrix:extract.ParameterSet]{*extract}()} returns.
}

\author{Henrik Bengtsson \url{http://www.braju.com/R/}}

\seealso{
  \code{\link[aroma.affymetrix:extract.ParameterSet]{*extract}()}.
  \code{\link[aroma.affymetrix:000-LB--LB-.ParameterSet]{*[[}()}.
  For more information see \code{\link{ParameterSet}}.
}
\keyword{internal}
\keyword{methods}
