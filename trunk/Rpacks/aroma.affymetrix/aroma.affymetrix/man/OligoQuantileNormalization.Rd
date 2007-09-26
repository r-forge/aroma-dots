%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  OligoQuantileNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{OligoQuantileNormalization}
\docType{class}
\alias{OligoQuantileNormalization}


 \title{The OligoQuantileNormalization class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class OligoQuantileNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{QuantileNormalization}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{OligoQuantileNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{OligoQuantileNormalization}\cr
extends \link[aroma.affymetrix]{QuantileNormalization}\cr



  This class represents a special \code{\link{QuantileNormalization}} that produces
  identical output as the quantile normalization in the \pkg{oligo} package.
 }

 \usage{OligoQuantileNormalization(...)}

 \arguments{
   \item{...}{Arguments passed to the constructor of
       \code{\link{QuantileNormalization}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from QuantileNormalization}:\cr
clearCache, process

 \bold{Methods inherited from Transform}:\cr
getFullName, getInputDataSet, getName, getOutputDataSet, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
