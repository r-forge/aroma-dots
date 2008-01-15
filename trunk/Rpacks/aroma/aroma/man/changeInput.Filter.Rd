%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Filter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{changeInput.Filter}
\alias{changeInput.Filter}
\alias{Filter.changeInput}
\alias{changeInput.Filter}
\alias{changeInput,Filter-method}


\title{Change input(s) on this filter and optionally all filters down the stream}

\usage{\method{changeInput}{Filter}(...)}

\description{
  Changes some input(s) to new input(s). A filter has one or more
  input connections that each can be connected to for instance  another
  filter's output or to some data objects. After having designed a
  network of filters connected to some data objects it is sometime
  desireble to connect the same network of filters to another set of
  data objects. This can be done easily by calling this method on the
  very last filter and tell it to recursively update all other filters
  accordingly.
}

\arguments{
  \item{newInput}{All connections to this input will be
    disconnected and replaced by the input object specified by
    this argument.}
  \item{oldInput}{All connections matching to this input will be replaced.
    If \code{\link[base]{NULL}} all connections will be replaced.}
  \item{recursive}{If \code{\link[base:logical]{TRUE}} this filter and all filters up the
    stream will update its connections according to \code{oldInput} and
   \code{newInput}. If \code{\link[base:logical]{FALSE}}, only this filter will update its
   input connection(s).}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  To the current set of connected inputs see \code{\link[aroma:getInput.Filter]{*getInput}()}.
  For more information see \code{\link{Filter}}..
}
\keyword{internal}
\keyword{methods}
