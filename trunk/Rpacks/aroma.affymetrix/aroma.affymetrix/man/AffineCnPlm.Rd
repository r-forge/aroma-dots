%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffineCnPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffineCnPlm}
\docType{class}
\alias{AffineCnPlm}


 \title{The AffineCnPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffineCnPlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffinePlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffineSnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{CnPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{AffineCnPlm}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AffineCnPlm}\cr
extends \link[aroma.affymetrix]{CnPlm}\cr


 }

 \usage{AffineCnPlm(..., combineAlleles=FALSE, tags="*")}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffineSnpPlm}}.}
   \item{combineAlleles}{If \code{\link[base:logical]{FALSE}}, allele A and allele B are treated
      seperately, otherwise together.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from CnPlm}:\cr
getCellIndices, getChipEffects, getProbeAffinities, setCombineAlleles

 \bold{Methods inherited from SnpPlm}:\cr
getCellIndices, getChipEffects, getProbeAffinities, setMergeStrands

 \bold{Methods inherited from AffinePlm}:\cr
getProbeAffinities

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
