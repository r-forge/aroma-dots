%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  OrFilter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{OrFilter}
\docType{class}
\alias{OrFilter}


\title{The OrFilter class}

\description{
 Package:  aroma \cr
\bold{Class OrFilter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{Filter}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma]{ParallelFilter}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{OrFilter}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{OrFilter}\cr
extends \emph{\link[aroma]{ParallelFilter}}\cr



  An OrFilter is a ParallelFilter that passes through indices that are
  outputs from \emph{any} input filters connected to this filter.
  This corresponds to the logical operator \emph{OR}. The filter could
  also have been called an \emph{union filter}.
}

\usage{OrFilter(...)}

\arguments{
  \item{...}{The input \code{\link{Filter}}s to be connected to.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma:getIndex.OrFilter]{getIndex}} \tab Gets indices accepted by this filter.\cr
}


\bold{Methods inherited from ParallelFilter}:\cr
changeInput, getInput

\bold{Methods inherited from Filter}:\cr
as.character, changeInput, getIndex, getInput, getParameter, getVisible, highlight, setParameter, setVisible, text

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

}


\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{\dontrun{See help(ParallelFilter) for an example.}}

\seealso{
  \code{\link{AndFilter}} and \code{\link{ParallelFilter}}.
}
\keyword{classes}
