%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% on Sat Jan 06 01:47:52 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{lapply.AffymetrixFileSet}
\alias{lapply.AffymetrixFileSet}
\alias{AffymetrixFileSet.lapply}
\alias{lapply.AffymetrixFileSet}
\alias{lapply,AffymetrixFileSet-method}


\title{Applies a function to each of the Affymetrix file files}

\description{
  Applies a function to each of the Affymetrix file files.
}

\usage{\method{lapply}{AffymetrixFileSet}(this, ...)}

\arguments{
 \item{...}{Not used.}
}

\value{
  Returns a named \code{\link[base]{list}} where the names are the names of the date files
  according to \code{\link[aroma.affymetrix:getNames.AffymetrixFileSet]{*getNames}()}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixFileSet}}.
}
\keyword{internal}
\keyword{methods}
