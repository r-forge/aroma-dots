%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PcrProbeCalibration.R
% 
% on Sat Jan 06 01:54:34 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{PcrProbeCalibration}
\docType{class}
\alias{PcrProbeCalibration}


 \title{The PcrProbeCalibration class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class PcrProbeCalibration}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Preprocessing}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbePreprocessing}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{PcrProbeCalibration}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{PcrProbeCalibration}\cr
extends \emph{\link[aroma.affymetrix]{ProbePreprocessing}}\cr



  This class represents the first step of pre-processing that corrects for
  PCR effects in the \emph{probe-level} intensities due to probe sequences
  and PCR fragment lengths, cf. Carvalho et al (2006).
 }

 \usage{PcrProbeCalibration(..., subsetToFit=1/5, flavor=c("*", "oligo"))}

 \arguments{
   \item{...}{Arguments passed to the constructor of
      \code{\link{ProbePreprocessing}}.}
   \item{subsetToFit}{The fraction of PM probes to be used for fitting
     the calibration function.}
   \item{flavor}{A \code{\link[base]{character}} string specifying what method to immitate.
     If \code{"oligo"}, the calibration gives identical results to the
     \pkg{oligo} package; the \code{subsetToFit} argument is ignored.
   }
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{process} \tab  -\cr
}


\bold{Methods inherited from ProbePreprocessing}:\cr
getRootPath

\bold{Methods inherited from Preprocessing}:\cr
as.character, getAsteriskTag, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputFiles, getParameters, getParametersAsString, getPath, getRootPath, getTags, isDone, process, setTags

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getInternalAddress, getStaticInstance, hasField, hashCode, isReferable, ll, load, names, newInstance, novirtual, objectSize, print, save, staticCode


 }

 \details{
   For sample \eqn{i=1,...,I} and probe \eqn{k=1,...,K}, the PM intensity
   is modelled as:
   \deqn{
     \log_2(PM_{ik}) = \mu_{ik} + g(L_k) + \sum_{t=1}^{25}
                       \sum_{b} h_b(t) I(b_{kt}=b) + \xi_{ik},
   }
   where  \eqn{\xi} is \eqn{N(0,\sigma^2)}.
   The \eqn{g(L_k)} is the PCR fragment-length effect for a probe in a
   SNP of length \eqn{L_k}.  Thus, all probes in a SNP have the same
   fragment-length effect \eqn{g(L_k)} (independ on sample).
   The \eqn{h_b(t)} is the effect of having base \eqn{b=\{A,C,G,T\}}
   at nucleotide position \eqn{t=1,...,25}.
   The \eqn{\mu_{ik}} is the signal of interest.
 }

 \section{Relationship to the oligo package}{
  The methods of class are adopted from the \code{\link[oligo]{snprma}} function
  in the \pkg{oligo} package. By setting argument
  \code{subsetToFit="asOligo"}, this calibration method produces the same
  result as the \pkg{oligo} package.
 }

 \section{Requirments}{
   To use this class the platform-design package for the chip type of the
   data set as well as the \pkg{oligo} package must be is installed.
   The methods, which are adopted from the \pkg{oligo} package, have been
   rewritten to minimize the memory usage. For instance, the platform-design
   package is never loaded, but instead subsets of its information is queried
   on demand utilizing the \code{\link{PlatformDesign}}.
 }

 \examples{\dontrun{
 }}

 \author{
   Henrik Bengtsson.
   All credits to Benilton Carvalho and Rafael Irizarry for the implementation
   of the underlying calibration algorithm.
 }

 \references{
   Carvalho B, Bengtsson H, Speed TP and Irizarry RA. \emph{Exploration,
   Normalization, and Genotype Calls of High Density Oligonucleotide SNP
   Array Data}, 2006.\cr

   Naef F and Magnasco MO. \emph{Solving the riddle of the bright
   mismatches: labeling and effective binding in oligonucleotide arrays}.
   Phys Rev E Stat Nonlin Soft Matter Phys, 2003, 68.\cr
 }
\keyword{classes}
