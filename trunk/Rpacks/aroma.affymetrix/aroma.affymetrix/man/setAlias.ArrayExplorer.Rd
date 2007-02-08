%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ArrayExplorer.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setAlias.ArrayExplorer}
\alias{setAlias.ArrayExplorer}
\alias{ArrayExplorer.setAlias}
\alias{setAlias.ArrayExplorer}
\alias{setAlias,ArrayExplorer-method}


\title{Sets the alias of the output set}

\description{
  Sets the alias of the output set.
  If specified, the alias overrides the data set name, which is used by
  default.
}

\usage{\method{setAlias}{ArrayExplorer}(this, alias=NULL, ...)}

\arguments{
 \item{alias}{A \code{\link[base]{character}} string for the new alias of the output set.
  The alias must consists of valid filename characters, and must not
  contain commas, which are used to separate tags.}
 \item{...}{Not used.}
}

\value{
  Returns nothing.
}

\details{
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:getAlias.ArrayExplorer]{*getAlias}()}.
  \code{\link[aroma.affymetrix:getName.ArrayExplorer]{*getName}()}.
  For more information see \code{\link{ArrayExplorer}}.
}
\keyword{internal}
\keyword{methods}
