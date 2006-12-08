%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChipEffectPreprocessing.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ChipEffectPreprocessing}
\docType{class}
\alias{ChipEffectPreprocessing}


 \title{The ChipEffectPreprocessing class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ChipEffectPreprocessing}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Preprocessing}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{ChipEffectPreprocessing}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{FragmentLengthNormalization}\cr

public abstract static class \bold{ChipEffectPreprocessing}\cr
extends \emph{\link[aroma.affymetrix]{Preprocessing}}\cr



  This abstract class represents a preprocessor that transforms
  chip effects obtained from probe-level modelling.
 }

 \usage{ChipEffectPreprocessing(dataSet=NULL, ...)}

 \arguments{
   \item{dataSet}{The input data set as an \code{\link{ChipEffectSet}}.}
   \item{...}{Arguments passed to the constructor of \code{\link{Preprocessing}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getRootPath} \tab  -\cr
}


\bold{Methods inherited from Preprocessing}:\cr
as.character, getAsteriskTag, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputFiles, getParameters, getParametersAsString, getPath, getRootPath, getTags, isDone, process, setTags

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \details{
   Subclasses must implement the \code{process()} method.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
