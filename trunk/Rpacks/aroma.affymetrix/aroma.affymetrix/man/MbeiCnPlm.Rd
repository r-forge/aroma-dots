%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MbeiCnPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{MbeiCnPlm}
\docType{class}
\alias{MbeiCnPlm}


 \title{The MbeiCnPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class MbeiCnPlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{UnitGroupsModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{MbeiPlm}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{MbeiSnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{CnPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{MbeiCnPlm}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{MbeiCnPlm}\cr
extends \link[aroma.affymetrix]{CnPlm}\cr


 }

 \usage{MbeiCnPlm(..., combineAlleles=FALSE, name="modelMbeiCnPlm")}

 \arguments{
   \item{...}{Arguments passed to \code{\link{MbeiSnpPlm}}.}
   \item{combineAlleles}{If \code{\link[base:logical]{FALSE}}, allele A and allele B are treated
      seperately, otherwise together.}
   \item{name}{The name of the PLM, also used as part of the pathname.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


\bold{Methods inherited from CnPlm}:\cr
getCellIndices, getChipEffects, getChipEffectSetClass, getFitUnitFunction, getProbeAffinities, getSubname, setCombineAlleles

\bold{Methods inherited from SnpPlm}:\cr
getCellIndices, getChipEffects, getChipEffectSetClass, getProbeAffinities, getSubname, setMergeStrands

\bold{Methods inherited from MbeiSnpPlm}:\cr


\bold{Methods inherited from MbeiPlm}:\cr
getFitFunction, getProbeAffinities

\bold{Methods inherited from ProbeLevelModel}:\cr
findUnitsTodo, fit, getCellIndices, getChipEffects, getChipEffectSetClass, getFitFunction, getFitUnitFunction, getProbeAffinities, highlight, plotMvsPosition, readUnits

\bold{Methods inherited from UnitGroupsModel}:\cr
as.character, fit, getCdf, getDataSet, getLabel, getName, getParameterSet, getPath, setLabel, setup

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
