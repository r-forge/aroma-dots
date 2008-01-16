%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AndFilter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AndFilter}
\docType{class}
\alias{AndFilter}


\title{The AndFilter class}

\description{
 Package:  aroma \cr
\bold{Class AndFilter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{Filter}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma]{ParallelFilter}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{AndFilter}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AndFilter}\cr
extends \emph{\link[aroma]{ParallelFilter}}\cr



  An AndFilter is a ParallelFilter that passes through indices that are
  outputs from \emph{all} input filters connected to this filter.
  This corresponds to the logical operator \emph{AND}. The filter could
  also have been called an \emph{intersection filter}.
}

\usage{AndFilter(...)}

\arguments{
  \item{...}{The input \code{\link{Filter}}s to be connected to.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma:getIndex.AndFilter]{getIndex}} \tab Gets indices accepted by this filter.\cr
}


 \bold{Methods inherited from ParallelFilter}:\cr
changeInput, getInput

 \bold{Methods inherited from Filter}:\cr
as.character, changeInput, getIndex, getInput, getParameter, getVisible, highlight, setParameter, setVisible, text

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}


\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{\dontrun{See help(ParallelFilter) for an example.}}

\seealso{
  \code{\link{OrFilter}} and \code{\link{ParallelFilter}}.
}
\keyword{classes}
