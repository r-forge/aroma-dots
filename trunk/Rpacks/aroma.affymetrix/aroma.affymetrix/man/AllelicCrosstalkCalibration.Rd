%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AllelicCrosstalkCalibration.R
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

 \usage{AllelicCrosstalkCalibration(..., targetAvg=2200, alpha=c(0.1, 0.075, 0.05, 0.03, 0.01), q=2, Q=98)}

 \arguments{
   \item{...}{Arguments passed to the constructor of
     \code{\link{ProbeLevelTransform}}.}
   \item{targetAvg}{The signal that average allele A and average allele B
     signals should have after calibration.}
   \item{alpha, q, Q}{}
 }

 \section{What probe signals are updated}{
    Calibration for crosstalk between allele signals applies by definition
    only SNP units.
    Note that, non-calibrated signals will be saved in the output files.
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
