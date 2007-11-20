%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSetReporter.R
% 
% on Tue Nov 20 12:39:31 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setAlias.AffymetrixFileSetReporter}
\alias{setAlias.AffymetrixFileSetReporter}
\alias{AffymetrixFileSetReporter.setAlias}
\alias{setAlias.AffymetrixFileSetReporter}
\alias{setAlias,AffymetrixFileSetReporter-method}


\title{Sets the alias of the output set}

\description{
  Sets the alias of the output set.
  If specified, the alias overrides the data set name, which is used by
  default.
}

\usage{\method{setAlias}{AffymetrixFileSetReporter}(this, alias=NULL, ...)}

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
  \code{\link[aroma.affymetrix:getAlias.AffymetrixFileSetReporter]{*getAlias}()}.
  \code{\link[aroma.affymetrix:getName.AffymetrixFileSetReporter]{*getName}()}.
  For more information see \code{\link{AffymetrixFileSetReporter}}.
}
\keyword{internal}
\keyword{methods}
