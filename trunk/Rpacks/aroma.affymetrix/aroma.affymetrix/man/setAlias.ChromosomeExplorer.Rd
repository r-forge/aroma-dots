%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChromosomeExplorer.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setAlias.ChromosomeExplorer}
\alias{setAlias.ChromosomeExplorer}
\alias{ChromosomeExplorer.setAlias}
\alias{setAlias.ChromosomeExplorer}
\alias{setAlias,ChromosomeExplorer-method}


\title{Sets the alias of the output set}

\description{
  Sets the alias of the output set.
  If specified, the alias overrides the model name, which is used by
  default.
}

\usage{\method{setAlias}{ChromosomeExplorer}(this, alias=NULL, ...)}

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
  \code{\link[aroma.affymetrix:getAlias.ChromosomeExplorer]{*getAlias}()}.
  \code{\link[aroma.affymetrix:getName.ChromosomeExplorer]{*getName}()}.
  For more information see \code{\link{ChromosomeExplorer}}.
}
\keyword{internal}
\keyword{methods}
