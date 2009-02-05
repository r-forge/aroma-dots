%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RocData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{RocData}
\docType{class}
\alias{RocData}


 \title{The RocData class}

 \description{
  Package:  aroma.cn.eval \cr
\bold{Class RocData}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{RocData}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.cn.eval]{MatrixRocData}\cr

public static class \bold{RocData}\cr
extends \link[R.oo]{Object}\cr


 }

 \usage{RocData(truth=NULL, data=NULL, recall=NULL, ...)}

 \arguments{
   \item{truth}{A \code{\link[base]{numeric}}al \code{\link[base]{vector}} of length N.}
   \item{data}{A \code{\link[base]{numeric}}al object of length N.}
   \item{recall}{(Optional) Unless \code{truth} is given as binary
     \eqn{\{0,1\}} values, it can be reclassified as such.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{anyMissing} \tab  -\cr
 \tab \code{as.character} \tab  -\cr
 \tab \code{as.data.frame} \tab  -\cr
 \tab \code{clearCache} \tab  -\cr
 \tab \code{dim} \tab  -\cr
 \tab \code{extractDataFrame} \tab  -\cr
 \tab \code{findTpAtFp} \tab  -\cr
 \tab \code{fit} \tab  -\cr
 \tab \code{getCallRate} \tab  -\cr
 \tab \code{getChecksum} \tab  -\cr
 \tab \code{getData} \tab  -\cr
 \tab \code{getInternalRocFields} \tab  -\cr
 \tab \code{getOrder} \tab  -\cr
 \tab \code{getRocFields} \tab  -\cr
 \tab \code{getTruth} \tab  -\cr
 \tab \code{isMissing} \tab  -\cr
 \tab \code{loadCachedResult} \tab  -\cr
 \tab \code{nbrOfDataPoints} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \details{
   Results are memoized, that is, cached to file, in order to speed up
   subsequent processing.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
