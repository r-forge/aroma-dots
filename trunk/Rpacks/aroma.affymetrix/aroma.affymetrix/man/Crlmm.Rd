%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Crlmm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{Crlmm}
\docType{class}
\alias{Crlmm}


 \title{The Crlmm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class Crlmm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{Crlmm}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{Crlmm}\cr
extends \link[R.oo]{Object}\cr


 }

 \usage{Crlmm(ces=NULL, minLLRforCalls=c(AA = 50, AB = 40, BB = 50), transform=c("log", "asinh"), tags="CRLMM", ...)}

 \arguments{
   \item{ces}{A \code{\link{SnpChipEffectSet}} to which this model should
     be fitted.}
   \item{minLLRforCalls}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length three.}
   \item{transform}{The transform used for the chip effects.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be appended to the tags of
      the input data set.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fit.Crlmm]{fit}} \tab Fits the model.\cr
 \tab \code{\link[aroma.affymetrix:getCdf.Crlmm]{getCdf}} \tab Gets the CDF structure for this model.\cr
 \tab \code{\link[aroma.affymetrix:getChipEffects.Crlmm]{getChipEffects}} \tab Gets the chip effects for this model.\cr
 \tab \code{\link[aroma.affymetrix:getFullName.Crlmm]{getFullName}} \tab Gets the full name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getName.Crlmm]{getName}} \tab Gets the name of the output summarized data set.\cr
 \tab \code{\link[aroma.affymetrix:getPath.Crlmm]{getPath}} \tab Gets the path of this model.\cr
 \tab \code{\link[aroma.affymetrix:getRootPath.Crlmm]{getRootPath}} \tab Gets the root path of this model.\cr
 \tab \code{\link[aroma.affymetrix:getTags.Crlmm]{getTags}} \tab Gets the tags of the output data set.\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
 }
\keyword{classes}
