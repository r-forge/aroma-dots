%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Filter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getParameter.Filter}
\alias{getParameter.Filter}
\alias{Filter.getParameter}
\alias{getParameter.Filter}
\alias{getParameter,Filter-method}


\title{Gets the values for a specific parameter for indices accepted by the filter}

\usage{\method{getParameter}{Filter}(this, paramName, ...)}

\description{
  Gets the values for a specific parameter for those indices that the
  filter lets through. If the filter's has this parameter set, e.g.
  through \code{\link[aroma:setParameter.Filter]{*setParameter}()}, all indices accepted by the filter will
  be having this parameter set to the filters value.
  If the parameter is not set in the filter, or its value is \code{\link[base]{NULL}}
  or \\code{\link[base]{NA}}, the filter asks all its input connections for assign
  parameters to the indices. If more than one input filter tries to set
  the parameter of an index, only the \emph{last} assignment will be used.
  Hence, the order in which the input filters are specified might affect
  the final assignment of paramters.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  To set a parameter of a filter see \code{\link[aroma:setParameter.Filter]{*setParameter}()}.
  To get the indices accepted by a filter see \code{\link[aroma:getIndex.Filter]{*getIndex}()}.
  For more information see \code{\link{Filter}}..
}
\keyword{internal}
\keyword{methods}
