%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Filter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{highlight.Filter}
\alias{highlight.Filter}
\alias{Filter.highlight}
\alias{highlight.Filter}
\alias{highlight,Filter-method}


\title{Highlights the data points accepted by the filter}

\usage{\method{highlight}{Filter}(this, cex="filter", col="filter", pch="filter", ..., recursive=FALSE)}

\description{
  Highlights the data points accepted by the filter and optionally its
  input objects.
  If the argument \code{recursive == TRUE} the filter will \emph{first}
  ask its input filters to highlight the data accepted by them. If these
  filters have filters connected to them they will before doing anything
  asking there input filters and so on. The last filter to highlight the
  data will therefore be this filter.
}

\arguments{
  \item{cex}{The scaling factor for \emph{all} points, i.e. for
    \emph{all} indices. If \code{\link[base]{NULL}}, the global default value will be
    used. If \code{"filter"} the scaling factor set by the filter will
    be used.}
  \item{col}{The color for \emph{all} points, i.e. for
    \emph{all} indices. If \code{\link[base]{NULL}}, the global default value will be
    used. If \code{"filter"} the color set by the filter will
    be used.}
  \item{pch}{The point style for \emph{all} points, i.e. for
    \emph{all} indices. If \code{\link[base]{NULL}}, the global default value will be
    used. If \code{"filter"} the point style set by the filter will
    be used.}
  \item{...}{Any parameters accepted by most standard plot functions.}
  \item{recursive}{If \code{\link[base:logical]{TRUE}} this filter and all filters up the
    stream will highlight the data points that passes through them.
    If \code{\link[base:logical]{FALSE}}, only this filter will be used to highlight the data.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  To label data points see \code{\link[aroma:text.Filter]{*text}()}.
  For more information see \code{\link{Filter}}..
}
\keyword{internal}
\keyword{methods}
