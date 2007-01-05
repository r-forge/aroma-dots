%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbePreprocessing.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ProbePreprocessing}
\docType{class}
\alias{ProbePreprocessing}


 \title{The ProbePreprocessing class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ProbePreprocessing}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Preprocessing}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{ProbePreprocessing}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AllelicCrosstalkCalibration}, \link[aroma.affymetrix]{OligoQuantileNormalization}, \link[aroma.affymetrix]{PcrProbeCalibration}, \link[aroma.affymetrix]{QuantileNormalization}\cr

public abstract static class \bold{ProbePreprocessing}\cr
extends \emph{\link[aroma.affymetrix]{Preprocessing}}\cr



  This abstract class represents a preprocessor that transforms
  probe-level signals, typically intensities.
 }

 \usage{ProbePreprocessing(...)}

 \arguments{
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
