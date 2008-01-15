%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getExtra.MicroarrayData}
\alias{getExtra.MicroarrayData}
\alias{MicroarrayData.getExtra}
\alias{getExtra.MicroarrayData}
\alias{getExtra,MicroarrayData-method}


\title{Gets an extra field}

\usage{\method{getExtra}{MicroarrayData}(this, key, ...)}

\arguments{
  \item{key}{The name under which the extra field is stored.}
}

\description{
  Gets an optional and extra field.
}

\examples{
  raw <- RawData()
  setExtra(raw, "date", date())
  getExtra(raw, "date")
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
