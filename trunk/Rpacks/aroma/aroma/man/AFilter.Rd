%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AFilter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AFilter}
\docType{class}
\alias{AFilter}


\title{FieldFilter for the log-intensities (A)}

\description{
 Package:  aroma \cr
\bold{Class AFilter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{Filter}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma]{SerialFilter}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma]{FieldFilter}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{AFilter}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AFilter}\cr
extends \link[aroma]{FieldFilter}\cr



  FieldFilter for the log-intensities (A).
}

\usage{AFilter(input, ...)}

\arguments{
  \item{input}{The input \code{\link{MicroarrayData}} object.}
  \item{...}{Any arguments accepted by the \code{\link{FieldFilter}} constructor.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\emph{No methods defined}.


\bold{Methods inherited from FieldFilter}:\cr
getIndex, getIndex0

\bold{Methods inherited from SerialFilter}:\cr
changeInput, getInput

\bold{Methods inherited from Filter}:\cr
as.character, changeInput, getIndex, getInput, getParameter, getVisible, highlight, setParameter, setVisible, text

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

}


\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  See also the \code{\link{FieldFilter}} class.
}

\keyword{manip}
\keyword{classes}
