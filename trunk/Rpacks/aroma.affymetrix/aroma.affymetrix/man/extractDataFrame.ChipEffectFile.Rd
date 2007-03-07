%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectFile.extractDataFrame.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{extractDataFrame.ChipEffectFile}
\alias{extractDataFrame.ChipEffectFile}
\alias{ChipEffectFile.extractDataFrame}
\alias{extractDataFrame.ChipEffectFile}
\alias{extractDataFrame,ChipEffectFile-method}


\title{Extract data as a data frame}

\description{
 Extract data as a data frame.
}

\usage{\method{extractDataFrame}{ChipEffectFile}(this, ..., units=NULL, field=c("theta", "sdTheta"), verbose=FALSE)}

\arguments{
  \item{units}{(The subset of units to be matched.
    If \code{\link[base]{NULL}}, all units are considered.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an Jx4 \code{\link[base]{data.frame}} where J is the number of units.
 The names of the last column is the name of the array.
 The names of the rows are the unit indices (as indexed by the CDF).
 The rows are ordered according to \code{units} arguments.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ChipEffectFile}}.
}
\keyword{internal}
\keyword{methods}
