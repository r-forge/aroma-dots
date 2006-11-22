%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AllelicCrosstalkCalibrator.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AllelicCrosstalkCalibrator}
\docType{class}
\alias{AllelicCrosstalkCalibrator}


 \title{The AllelicCrosstalkCalibrator class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AllelicCrosstalkCalibrator}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{AllelicCrosstalkCalibrator}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AllelicCrosstalkCalibrator}\cr
extends \link[R.oo]{Object}\cr



  This class represents a calibration function that transforms the
  probe-level signals such that the signals from the two alleles are
  orthogonal.
 }

 \usage{AllelicCrosstalkCalibrator(dataSet=NULL, tags="*", ...)}

 \arguments{
   \item{dataSet}{A \code{\link{AffymetrixCelSet}}.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be appended to the tags of
      the input data set.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getFullName.AllelicCrosstalkCalibrator]{getFullName}} \tab Gets the full name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getInputDataSet.AllelicCrosstalkCalibrator]{getInputDataSet}} \tab Gets the source data set.\cr
 \tab \code{\link[aroma.affymetrix:getName.AllelicCrosstalkCalibrator]{getName}} \tab Gets the name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getOutputDataSet.AllelicCrosstalkCalibrator]{getOutputDataSet}} \tab Gets the normalized data set.\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getPath.AllelicCrosstalkCalibrator]{getPath}} \tab Gets the path of the output data set.\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getTags.AllelicCrosstalkCalibrator]{getTags}} \tab Gets the tags of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:isDone.AllelicCrosstalkCalibrator]{isDone}} \tab Checks if the data set is calibrated or not.\cr
 \tab \code{\link[aroma.affymetrix:process.AllelicCrosstalkCalibrator]{process}} \tab Calibrates the data set.\cr
 \tab \code{setTags} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
