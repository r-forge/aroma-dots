%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setExtra.MicroarrayData}
\alias{setExtra.MicroarrayData}
\alias{MicroarrayData.setExtra}
\alias{setExtra.MicroarrayData}
\alias{setExtra,MicroarrayData-method}


\title{Sets an extra field}

\usage{\method{setExtra}{MicroarrayData}(this, key, value, ...)}

\arguments{
  \item{key}{The name under which the extra field should be stored.}
  \item{value}{The value to be stored.}
}

\description{
  Stores optional and extra information of any format.
}

\examples{
  raw <- RawData()
  setExtra(raw, "date", date())
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
