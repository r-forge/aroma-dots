%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RmaCnPlm.R
% 
% on Sat Jan 06 01:56:29 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{RmaCnPlm}
\docType{class}
\alias{RmaCnPlm}


 \title{The RmaCnPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class RmaCnPlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitGroupsModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{RmaPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{RmaSnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{CnPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{RmaCnPlm}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{RmaCnPlm}\cr
extends \link[aroma.affymetrix]{CnPlm}\cr


 }

 \usage{RmaCnPlm(..., combineAlleles=FALSE, tags="*")}

 \arguments{
   \item{...}{Arguments passed to \code{\link{RmaSnpPlm}}.}
   \item{combineAlleles}{If \code{\link[base:logical]{FALSE}}, allele A and allele B are treated
      seperately, otherwise together.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


\bold{Methods inherited from CnPlm}:\cr
getCellIndices, getChipEffects, getChipEffectSetClass, getFitUnitFunction, getProbeAffinities, setCombineAlleles

\bold{Methods inherited from SnpPlm}:\cr
getCellIndices, getChipEffects, getChipEffectSetClass, getProbeAffinities, setMergeStrands

\bold{Methods inherited from RmaSnpPlm}:\cr


\bold{Methods inherited from RmaPlm}:\cr
getFitFunction, getProbeAffinities

\bold{Methods inherited from ProbeLevelModel}:\cr
findUnitsTodo, fit, getCellIndices, getChipEffects, getChipEffectSetClass, getFitFunction, getFitUnitFunction, getProbeAffinities, getRootPath, highlight, plotMvsPosition, readUnits

\bold{Methods inherited from UnitGroupsModel}:\cr
getLabel, setLabel, setup

\bold{Methods inherited from UnitModel}:\cr
as.character, findUnitsTodo, fit, getCdf, getCellIndices, getDataSet, getFitUnitFunction, getFullName, getName, getParameterSet, getPath, getRootPath, getTags, readUnits, setTags

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getInternalAddress, getStaticInstance, hasField, hashCode, isReferable, ll, load, names, newInstance, novirtual, objectSize, print, save, staticCode


 }

 \section{Model}{
   TO DO.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
