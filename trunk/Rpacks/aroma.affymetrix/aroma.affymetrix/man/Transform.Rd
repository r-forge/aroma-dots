%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Transform.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{Transform}
\docType{class}
\alias{Transform}


 \title{The Transform class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class Transform}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{Transform}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AllelicCrosstalkCalibration}, \emph{\link[aroma.affymetrix]{BackgroundCorrection}}, \link[aroma.affymetrix]{ChipEffectGroupMerge}, \emph{\link[aroma.affymetrix]{ChipEffectTransform}}, \link[aroma.affymetrix]{DChipQuantileNormalization}, \link[aroma.affymetrix]{FragmentEquivalentClassNormalization}, \link[aroma.affymetrix]{FragmentLengthNormalization}, \link[aroma.affymetrix]{GcContentNormalization}, \link[aroma.affymetrix]{GcRmaBackgroundCorrection}, \link[aroma.affymetrix]{OpticalBackgroundCorrection}, \emph{\link[aroma.affymetrix]{ProbeLevelTransform}}, \link[aroma.affymetrix]{QuantileNormalization}, \link[aroma.affymetrix]{RmaBackgroundCorrection}, \link[aroma.affymetrix]{ScaleNormalization}, \link[aroma.affymetrix]{SnpChipEffectGroupMerge}\cr

public abstract static class \bold{Transform}\cr
extends \link[R.oo]{Object}\cr



  This abstract class represents a transform (algorithm/operator) that
  transforms data.  A transform has an input data set, which is
  transformed into an output data set.
 }

 \usage{Transform(dataSet=NULL, tags="*", ...)}

 \arguments{
   \item{dataSet}{The input data set as an \code{\link{AffymetrixCelSet}}.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be appended to the tags of
      the input data set.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:getFullName.Transform]{getFullName}} \tab Gets the full name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getInputDataSet.Transform]{getInputDataSet}} \tab Gets the source data set.\cr
 \tab \code{\link[aroma.affymetrix:getName.Transform]{getName}} \tab Gets the name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getOutputDataSet.Transform]{getOutputDataSet}} \tab Gets the transformed data set.\cr
 \tab \code{\link[aroma.affymetrix:getPath.Transform]{getPath}} \tab Gets the path of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getTags.Transform]{getTags}} \tab Gets the tags of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:isDone.Transform]{isDone}} \tab Checks if the data set is processed or not.\cr
 \tab \code{\link[aroma.affymetrix:process.Transform]{process}} \tab Processes the data set.\cr
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
