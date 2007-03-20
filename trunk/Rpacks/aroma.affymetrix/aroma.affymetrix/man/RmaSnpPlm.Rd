%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RmaSnpPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{RmaSnpPlm}
\docType{class}
\alias{RmaSnpPlm}


 \title{The RmaSnpPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class RmaSnpPlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{RmaPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{RmaSnpPlm}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{RmaCnPlm}\cr

public static class \bold{RmaSnpPlm}\cr
extends \link[aroma.affymetrix]{SnpPlm}\cr


 }

 \usage{RmaSnpPlm(..., mergeStrands=FALSE, tags="*")}

 \arguments{
   \item{...}{Arguments passed to \code{\link{RmaPlm}}.}
   \item{mergeStrands}{If \code{\link[base:logical]{TRUE}}, the sense and the anti-sense strands are
      fitted together, otherwise separately.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from SnpPlm}:\cr
getCellIndices, getChipEffects, getProbeAffinities, setMergeStrands

 \bold{Methods inherited from RmaPlm}:\cr
getCalculateResidualsFunction

 \bold{Methods inherited from ProbeLevelModel}:\cr
calculateResiduals, calculateResidualSet, calculateWeights, fit, getCalculateResidualsFunction, getChipEffects, getChipEffectSet, getProbeAffinities, getProbeAffinityFile, getResidualSet, getWeightsSet

 \bold{Methods inherited from UnitModel}:\cr
findUnitsTodo, getCellIndices, readUnits

 \bold{Methods inherited from Model}:\cr
fit, getAlias, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
