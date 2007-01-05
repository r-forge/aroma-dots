%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% on Sat Jan 06 01:47:53 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{seq.AffymetrixFileSet}
\alias{seq.AffymetrixFileSet}
\alias{AffymetrixFileSet.seq}
\alias{seq.AffymetrixFileSet}
\alias{seq,AffymetrixFileSet-method}


\title{Gets an vector of file indices}

\description{
  Gets an vector of file indices.
}

\usage{\method{seq}{AffymetrixFileSet}(this, ...)}

\arguments{
 \item{...}{Not used.}
}

\value{
  Returns an \code{\link[base]{integer}} \code{\link[base]{vector}} in [1,N] where N is the number of files,
  or an empty vector if the file set is empty.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixFileSet}}.
}
\keyword{internal}
\keyword{methods}
