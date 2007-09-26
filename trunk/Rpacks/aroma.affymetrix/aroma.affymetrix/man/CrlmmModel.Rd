%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CrlmmModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CrlmmModel}
\docType{class}
\alias{CrlmmModel}


 \title{The CrlmmModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CrlmmModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{CrlmmModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CrlmmModel}\cr
extends \emph{\link[aroma.affymetrix]{Model}}\cr



  This class represents the CRLMM model of Carvalho et al. (2006).
 }

 \usage{CrlmmModel(ces=NULL, minLLRforCalls=c(AA = 50, AB = 40, BB = 50), transform=c("log", "asinh"), tags="CRLMM", ...)}

 \arguments{
   \item{ces}{A \code{\link{SnpChipEffectSet}} to which this model should
     be fitted.}
   \item{minLLRforCalls}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length three.}
   \item{transform}{The transform used for the chip effects.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be appended to the tags of
      the input data set.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fit} \tab  -\cr
 \tab \code{getChipEffects} \tab  -\cr
}


 \bold{Methods inherited from Model}:\cr
fit, getAlias, getAsteriskTag, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Requirments}{
 
   This implementation is partly adopted from and partly makes use of
   the implementation available in the \pkg{oligo} package.
   It was implemented so it can give identical results to the \pkg{oligo}
   implementation.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
  Carvalho B, Bengtsson H, Speed TP, and Irizarry RA. \emph{Exploration,
  Normalization, and Genotype Calls of High Density Oligonucleotide SNP
  Array Data}, Biostatistics, 2006.\cr
 }
\keyword{classes}
