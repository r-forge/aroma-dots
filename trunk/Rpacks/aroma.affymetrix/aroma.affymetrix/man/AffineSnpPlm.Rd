%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffineSnpPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffineSnpPlm}
\docType{class}
\alias{AffineSnpPlm}


 \title{The AffineSnpPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffineSnpPlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{UnitGroupsModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffinePlm}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpPlm}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{AffineSnpPlm}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{AffineSnpPlm}\cr
extends \link[aroma.affymetrix]{SnpPlm}\cr



  This class represents the log-additive model used in RMA.
  It can be used to fit the model on a \code{\link{AffymetrixCelSet}}.
 }

 \usage{AffineSnpPlm(..., name="modelAffineSnpPlm", mergeStrands=FALSE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffinePlm}}.}
   \item{name}{The name of the model, which is also used in the pathname.}
   \item{mergeStrands}{If \code{\link[base:logical]{TRUE}}, the sense and the anti-sense strands are
      fitted together, otherwise separately.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


\bold{Methods inherited from SnpPlm}:\cr
getCellIndices, getChipEffects, getChipEffectSetClass, getProbeAffinities, setMergeStrands

\bold{Methods inherited from AffinePlm}:\cr
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
