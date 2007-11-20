%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AllelicCrosstalkCalibration.R
% 
% on Tue Nov 20 12:39:36 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AllelicCrosstalkCalibration}
\docType{class}
\alias{AllelicCrosstalkCalibration}


 \title{The AllelicCrosstalkCalibration class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AllelicCrosstalkCalibration}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{AllelicCrosstalkCalibration}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AllelicCrosstalkCalibration}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelTransform}}\cr



  This class represents a calibration function that transforms the
  probe-level signals such that the signals from the two alleles are
  orthogonal.
  The method fits and calibrates PM signals only.  MM signals will not
  affect the model fitting and are unaffected.
 }

 \usage{AllelicCrosstalkCalibration(..., targetAvg=c(2200, 2200), subsetToAvg=NULL, alpha=c(0.1, 0.075, 0.05, 0.03, 0.01), q=2, Q=98)}

 \arguments{
   \item{...}{Arguments passed to the constructor of
     \code{\link{ProbeLevelTransform}}.}
   \item{targetAvg}{The signal(s) that either the average of the sum
     (if one target value) or the average of each of the alleles
     (if two target values) should have after calibration.}
   \item{subsetToAvg}{The indices of the cells (taken as the intersect of
     existing indices) used to calculate average in order to rescale to
     the target average. If \code{\link[base]{NULL}}, all probes are considered.}
   \item{alpha, q, Q}{}
 }

 \section{What probe signals are updated?}{
   Calibration for crosstalk between allele signals applies by definition
   only SNP units.  It is only PM probes that will be calibrated.
   Note that, non-calibrated signals will be saved in the output files.
 }

 \section{What probe signals are used to fit model?}{
   All PM probe pairs are used to fit the crosstalk model.
   In the second step where signals are rescaled to a target average,
   it is possible to specify the set of cells that should be included
   when estimating the target average.
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.AllelicCrosstalkCalibration]{process}} \tab Calibrates the data set.\cr
}


 \bold{Methods inherited from Transform}:\cr
getFullName, getInputDataSet, getName, getOutputDataSet, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
