%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Filter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{text.Filter}
\alias{text.Filter}
\alias{Filter.text}
\alias{text.Filter}
\alias{text,Filter-method}


\title{Labels the data points accepted by the filter}

\usage{\method{text}{Filter}(x, labels=NULL, cex="filter", col=NULL, ..., recursive=FALSE)}

\description{
  Labels the data points accepted by the filter and optionally its
  input objects.
  If the argument \code{recursive == TRUE} the filter will \emph{first}
  ask its input filters to label the data accepted by them. If these
  filters have filters connected to them they will before doing anything
  asking there input filters and so on. The last filter to label the
  data will therefore be this filter.
}

\arguments{
  \item{labels}{The labels for \emph{all} points, i.e. for \emph{all}
    indices. If \code{\link[base]{NULL}}, the indices will be used as labels.}
  \item{cex}{The scaling factor for \emph{all} points, i.e. for
    \emph{all} indices. If \code{\link[base]{NULL}}, the global default value will be
    used. If \code{"filter"} the scaling factor set by the filter will
    be used.}
  \item{col}{The color for \emph{all} points, i.e. for
    \emph{all} indices. If \code{\link[base]{NULL}}, the global default value will be
    used. If \code{"filter"} the color set by the filter will
    be used.}
  \item{...}{Any parameters accepted by standard \code{\link[base]{text}}.}
  \item{recursive}{If \code{\link[base:logical]{TRUE}} this filter and all filters up the
    stream will label the data points that passes through them.
    If \code{\link[base:logical]{FALSE}}, only this filter will be used to label the data.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  To highlight data points see \code{\link[aroma:highlight.Filter]{*highlight}()}.
  For more information see \code{\link{Filter}}..
}
\keyword{internal}
\keyword{methods}
