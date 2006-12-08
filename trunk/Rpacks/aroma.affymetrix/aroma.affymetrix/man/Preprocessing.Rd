%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Preprocessing.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{Preprocessing}
\docType{class}
\alias{Preprocessing}


 \title{The Preprocessing class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class Preprocessing}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{Preprocessing}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AllelicCrosstalkCalibration}, \emph{\link[aroma.affymetrix]{ChipEffectPreprocessing}}, \link[aroma.affymetrix]{FragmentLengthNormalization}, \link[aroma.affymetrix]{PcrProbeCalibration}, \emph{\link[aroma.affymetrix]{ProbePreprocessing}}, \link[aroma.affymetrix]{QuantileNormalization}\cr

public abstract static class \bold{Preprocessing}\cr
extends \link[R.oo]{Object}\cr



  This abstract class represents a preprocessor that fits a transformation
  function and transform the data set accordingly.
  A preprocessor has an input data set, which is transformed into an
  output data set.
 }

 \usage{Preprocessing(dataSet=NULL, tags="*", ...)}

 \arguments{
   \item{dataSet}{The input data set as an \code{\link{AffymetrixCelSet}}.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be appended to the tags of
      the input data set.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getFullName.Preprocessing]{getFullName}} \tab Gets the full name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getInputDataSet.Preprocessing]{getInputDataSet}} \tab Gets the source data set.\cr
 \tab \code{\link[aroma.affymetrix:getName.Preprocessing]{getName}} \tab Gets the name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getOutputDataSet.Preprocessing]{getOutputDataSet}} \tab Gets the preprocessed data set.\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getPath.Preprocessing]{getPath}} \tab Gets the path of the output data set.\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getTags.Preprocessing]{getTags}} \tab Gets the tags of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:isDone.Preprocessing]{isDone}} \tab Checks if the data set is processed or not.\cr
 \tab \code{\link[aroma.affymetrix:process.Preprocessing]{process}} \tab Processes the data set.\cr
 \tab \code{setTags} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \details{
   Subclasses must implement the \code{process()} method.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
