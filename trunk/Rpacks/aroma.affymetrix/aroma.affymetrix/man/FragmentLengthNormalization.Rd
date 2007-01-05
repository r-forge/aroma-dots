%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FragmentLengthNormalization.R
% 
% on Sat Jan 06 01:51:29 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{FragmentLengthNormalization}
\docType{class}
\alias{FragmentLengthNormalization}


 \title{The FragmentLengthNormalization class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class FragmentLengthNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Preprocessing}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectPreprocessing}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{FragmentLengthNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{FragmentLengthNormalization}\cr
extends \emph{\link[aroma.affymetrix]{ChipEffectPreprocessing}}\cr



  This class represents a normalization method that corrects for PCR
  fragment length effects on total copy-number chip-effect estimates.
 }

 \usage{FragmentLengthNormalization(dataSet=NULL, ..., targetFunction=NULL, subsetToFit=NULL)}

 \arguments{
   \item{dataSet}{A \code{\link{CnChipEffectSet}}.}
   \item{...}{Additional arguments passed to the constructor of
     \code{\link{ChipEffectPreprocessing}}.}
   \item{targetFunction}{A \code{\link[base]{function}}.  The target function to which all arrays
     should be normalized to.}
   \item{subsetToFit}{The units from which the normalization curve should
     be estimated.  If \code{\link[base]{NULL}}, all are considered.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCdf} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:process.FragmentLengthNormalization]{process}} \tab Normalizes the data set.\cr
}


\bold{Methods inherited from ChipEffectPreprocessing}:\cr
getRootPath

\bold{Methods inherited from Preprocessing}:\cr
as.character, getAsteriskTag, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputFiles, getParameters, getParametersAsString, getPath, getRootPath, getTags, isDone, process, setTags

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getInternalAddress, getStaticInstance, hasField, hashCode, isReferable, ll, load, names, newInstance, novirtual, objectSize, print, save, staticCode


 }

 \examples{\dontrun{
 }}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
