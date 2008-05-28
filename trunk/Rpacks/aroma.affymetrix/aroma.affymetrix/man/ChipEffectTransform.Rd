%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectTransform.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ChipEffectTransform}
\docType{class}
\alias{ChipEffectTransform}


 \title{The ChipEffectTransform class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ChipEffectTransform}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{ChipEffectTransform}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{ChipEffectGroupMerge}, \link[aroma.affymetrix]{FragmentEquivalentClassNormalization}, \link[aroma.affymetrix]{FragmentLengthNormalization}, \link[aroma.affymetrix]{GcContentNormalization}, \link[aroma.affymetrix]{SnpChipEffectGroupMerge}\cr

public abstract static class \bold{ChipEffectTransform}\cr
extends \emph{\link[aroma.affymetrix]{Transform}}\cr



  This abstract class represents a transform that transforms chip-effect
  estimates obtained from probe-level modelling.
 }

 \usage{ChipEffectTransform(dataSet=NULL, ...)}

 \arguments{
   \item{dataSet}{The input data set as an \code{\link{ChipEffectSet}}.}
   \item{...}{Arguments passed to the constructor of \code{\link{Transform}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getOutputDataSet} \tab  -\cr
}


 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getFullName, getInputDataSet, getName, getOutputDataSet, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \details{
   Subclasses must implement the \code{process()} method.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
