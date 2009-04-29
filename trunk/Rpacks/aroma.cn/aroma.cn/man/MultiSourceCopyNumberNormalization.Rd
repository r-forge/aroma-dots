%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MultiSourceCopyNumberNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{MultiSourceCopyNumberNormalization}
\docType{class}
\alias{MultiSourceCopyNumberNormalization}


 \title{The MultiSourceCopyNumberNormalization class}

 \description{
  Package:  aroma.cn \cr
\bold{Class MultiSourceCopyNumberNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{MultiSourceCopyNumberNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{MultiSourceCopyNumberNormalization}\cr
extends \link[R.oo]{Object}\cr



  A normalization method that normalizes copy-number estimates measured
  by multiple sites and/or platforms for common samples.  It normalizes the
  estimates toward a common scale such that for any copy-number level
  the mean level of the normalized data are the same.
 }

 \usage{MultiSourceCopyNumberNormalization(dsList=NULL, fitUgp=NULL, subsetToFit=NULL, targetDimension=1, ...)}

 \arguments{
  \item{dsList}{A \code{\link[base]{list}} of K \code{\link[aroma.core]{AromaUnitTotalCnBinarySet}}:s.}
  \item{fitUgp}{An \code{\link[aroma.core]{AromaUgpFile}} that specifies the
    common set of loci used to normalize the data sets at.}
  \item{subsetToFit}{The subset of loci (as mapped by the \code{fitUgp}
    object) to be used to fit the normalization functions.
    If \code{\link[base]{NULL}}, loci on chromosomes 1-22 are used, but not on ChrX and ChrY.
  }
  \item{targetDimension}{A \code{\link[base]{numeric}} index specifying the data set in
    \code{dsList} to which each platform in standardize towards.
    If \code{\link[base]{NULL}}, the arbitrary scale along the fitted principal curve
    is used.  This always starts at zero and increases.}
  \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.cn:getAllNames.MultiSourceCopyNumberNormalization]{getAllNames}} \tab Gets the names of all unique samples across all sources.\cr
 \tab \code{\link[aroma.cn:getInputDataSets.MultiSourceCopyNumberNormalization]{getInputDataSets}} \tab Gets the list of data sets to be normalized.\cr
 \tab \code{getOutputDataSets} \tab  -\cr
 \tab \code{getOutputPaths} \tab  -\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{nbrOfDataSets} \tab  -\cr
 \tab \code{\link[aroma.cn:process.MultiSourceCopyNumberNormalization]{process}} \tab Normalizes all samples.\cr
}


 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \details{
   The multi-source normalization method is by nature a single-sample method,
   that is, it normalizes arrays for one sample at the time and independently
   of all other samples/arrays.

   However, the current implementation is such that it first generates
   smoothed data for \emph{all} samples/arrays.  Then, it normalizes the
   sample one by one.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
