%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Filter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getVisible.Filter}
\alias{getVisible.Filter}
\alias{Filter.getVisible}
\alias{getVisible.Filter}
\alias{getVisible,Filter-method}


\title{Checks if the filter is visible or not}

\usage{\method{getVisible}{Filter}(obj, ...)}

\description{
  Checks if the filter is visible or not. A filter can be set to visible
  or non-visible if the data for the indices outputted by the filter
  should be highlighted, labeled etc. By default a filter is visible.
}

\value{
  Returns \code{\link[base:logical]{TRUE}} if the filter is visible, otherwise \code{\link[base:logical]{FALSE}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:setVisible.Filter]{*setVisible}()}.
  For more information see \code{\link{Filter}}..
}
\keyword{internal}
\keyword{methods}
