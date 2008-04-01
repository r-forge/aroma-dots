%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CbsModel.fitOne.R
% 
% on Mon Mar 31 17:07:02 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fitOne.CbsModel}
\alias{fitOne.CbsModel}
\alias{CbsModel.fitOne}
\alias{fitOne.CbsModel}
\alias{fitOne,CbsModel-method}


\title{Fits the CBS model for one chromosome in one sample}

\description{
 Fits the CBS model for one chromosome in one sample.
}

\usage{\method{fitOne}{CbsModel}(this, data, chromosome, ..., verbose=FALSE)}

\arguments{
  \item{data}{}
  \item{...}{}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the \code{\link[DNAcopy]{DNAcopy}} object returned by
 \code{\link[DNAcopy]{segment}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally \code{\link[DNAcopy]{segment}} is used.
  For more information see \code{\link{CbsModel}}.
}
\keyword{internal}
\keyword{methods}
