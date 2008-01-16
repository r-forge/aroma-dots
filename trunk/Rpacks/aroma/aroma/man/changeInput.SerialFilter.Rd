%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SerialFilter.R
% 
% on Tue Jan 15 18:37:04 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{changeInput.SerialFilter}
\alias{changeInput.SerialFilter}
\alias{SerialFilter.changeInput}
\alias{changeInput.SerialFilter}
\alias{changeInput,SerialFilter-method}


\title{Change input(s) on this filter and optionally all filters down the stream}

\usage{\method{changeInput}{SerialFilter}(this, newInput, oldInput=NULL, recursive=FALSE, ...)}

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
    If \code{\link[base]{NULL}} all connections will be replaced. Default value is
    \code{\link[base]{NULL}}.}
  \item{recursive}{If \code{\link[base:logical]{TRUE}} this filter and all filters up the
    stream will update its connections according to \code{oldInput} and
   \code{newInput}. If \code{\link[base:logical]{FALSE}}, only this filter will update its
   input connection(s).}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  To get the current set of connected inputs see \code{\link{Filter.getInput}}.
  For more information see \code{\link{SerialFilter}}.
}
\keyword{internal}
\keyword{methods}
