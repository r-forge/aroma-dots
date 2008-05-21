%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeLevelTransform.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ProbeLevelTransform}
\docType{class}
\alias{ProbeLevelTransform}


 \title{The ProbeLevelTransform class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ProbeLevelTransform}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{ProbeLevelTransform}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AllelicCrosstalkCalibration}, \emph{\link[aroma.affymetrix]{BackgroundCorrection}}, \link[aroma.affymetrix]{DChipQuantileNormalization}, \link[aroma.affymetrix]{GcRmaBackgroundCorrection}, \link[aroma.affymetrix]{OpticalBackgroundCorrection}, \link[aroma.affymetrix]{QuantileNormalization}, \link[aroma.affymetrix]{RmaBackgroundCorrection}, \link[aroma.affymetrix]{ScaleNormalization}, \link[aroma.affymetrix]{SpatialRowColumnNormalization}\cr

public abstract static class \bold{ProbeLevelTransform}\cr
extends \emph{\link[aroma.affymetrix]{Transform}}\cr



  This abstract class represents a transformation methods that transforms
  probe-level signals, typically intensities.
 }

 \usage{ProbeLevelTransform(...)}

 \arguments{
   \item{...}{Arguments passed to the constructor of \code{\link{Transform}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from Transform}:\cr
getFullName, getInputDataSet, getName, getOutputDataSet, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \details{
   Subclasses must implement the \code{process()} method.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
