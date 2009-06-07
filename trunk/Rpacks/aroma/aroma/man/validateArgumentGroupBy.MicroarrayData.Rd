%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.ARGS.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{validateArgumentGroupBy.MicroarrayData}
\alias{validateArgumentGroupBy.MicroarrayData}
\alias{MicroarrayData.validateArgumentGroupBy}
\alias{validateArgumentGroupBy.MicroarrayData}
\alias{validateArgumentGroupBy,MicroarrayData-method}


\title{Validates the argument 'groupBy'}

\description{
 Validates the argument 'groupBy'.
}

\usage{\method{validateArgumentGroupBy}{MicroarrayData}(this, groupBy, ...)}

\arguments{
  \item{groupBy}{Either a \code{\link{LayoutGroups}} object or a \code{\link[base]{character}}
    string specifying a layout group. If \code{\link[base]{NULL}}, all spots are
    considered.}
}

\value{
 Returns a \code{\link[base]{list}} containing elements that each contains a \code{\link[base]{vector}} of
 spot indices. If \code{\link[base]{NULL}}, the list will contain one element with
 \emph{all} indices.

 If the argument was invalid an \code{\link[R.oo]{Exception}} is thrown.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
