%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{seq.AffymetrixDataSet}
\alias{seq.AffymetrixDataSet}
\alias{AffymetrixDataSet.seq}
\alias{seq.AffymetrixDataSet}
\alias{seq,AffymetrixDataSet-method}


\title{Gets an vector of data file indices}

\description{
  Gets an vector of data file indices.
}

\usage{\method{seq}{AffymetrixDataSet}(this, ...)}

\arguments{
 \item{...}{Not used.}
}

\value{
  Returns an \code{\link[base]{integer}} \code{\link[base]{vector}} in [1,N] where N is the number of arrays,
  or an empty vector if the data set is empty.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixDataSet}}.
}
\keyword{internal}
\keyword{methods}
