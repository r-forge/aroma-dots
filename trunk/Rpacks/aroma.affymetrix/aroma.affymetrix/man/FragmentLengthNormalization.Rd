%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FragmentLengthNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{FragmentLengthNormalization}
\docType{class}
\alias{FragmentLengthNormalization}


 \title{The FragmentLengthNormalization class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class FragmentLengthNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{Transform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ChipEffectTransform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{FragmentLengthNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{FragmentLengthNormalization}\cr
extends \emph{\link[aroma.affymetrix]{ChipEffectTransform}}\cr



  This class represents a normalization method that corrects for PCR
  fragment length effects on copy-number chip-effect estimates.
 }

 \usage{FragmentLengthNormalization(dataSet=NULL, ..., targetFunctions=NULL, subsetToFit="-XY", shift=0)}

 \arguments{
   \item{dataSet}{A \code{\link{SnpChipEffectSet}}.}
   \item{...}{Additional arguments passed to the constructor of
     \code{\link{ChipEffectTransform}}.}
   \item{targetFunctions}{An optional list of \code{\link[base]{function}}s.
     For each enzyme there is one target function to which all arrays
     should be normalized to.}
   \item{subsetToFit}{The units from which the normalization curve should
     be estimated.  If \code{\link[base]{NULL}}, all are considered.}
   \item{shift}{An optional amount the data points should be shifted
      (translated).}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{clearCache} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:process.FragmentLengthNormalization]{process}} \tab Normalizes the data set.\cr
}


 \bold{Methods inherited from ChipEffectTransform}:\cr
getOutputDataSet

 \bold{Methods inherited from Transform}:\cr
getOutputDataSet, getOutputFiles

 \bold{Methods inherited from AromaTransform}:\cr
getFullName, getInputDataSet, getName, getOutputDataSet, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Requirements}{
   This class requires a SNP information annotation file for the
   chip type to be normalized.
 }

 \details{
   For SNPs, the normalization function is estimated based on the total
   chip effects, i.e. the sum of the allele signals.  The normalizing
   is done by rescale the chip effects on the intensity scale such that
   the mean of the total chip effects are the same across samples for
   any given fragment length.  For allele-specific estimates, both alleles
   are always rescaled by the same amount.  Thus, when normalizing
   allele-specific chip effects, the total chip effects is change, but not
   the relative allele signal, e.g. the allele B frequency.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
