%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setAlias.AffymetrixFileSet}
\alias{setAlias.AffymetrixFileSet}
\alias{AffymetrixFileSet.setAlias}
\alias{setAlias.AffymetrixFileSet}
\alias{setAlias,AffymetrixFileSet-method}


\title{Sets the alias of the file set}

\description{
  Sets the alias of the file set.
  If specified, the alias overrides the name inferred from the pathname
  of the file set.  This can be used in order to use another name of the
  output data set than the input data set of many transforms and models.
}

\usage{\method{setAlias}{AffymetrixFileSet}(this, alias=NULL, ...)}

\arguments{
 \item{alias}{A \code{\link[base]{character}} string for the new alias of the file set.
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
  \code{\link[aroma.affymetrix:getAlias.AffymetrixFileSet]{*getAlias}()}.
  \code{\link[aroma.affymetrix:getName.AffymetrixFileSet]{*getName}()}.
  For more information see \code{\link{AffymetrixFileSet}}.
}
\keyword{internal}
\keyword{methods}
