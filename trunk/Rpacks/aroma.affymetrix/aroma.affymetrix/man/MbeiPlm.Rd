%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MbeiPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{MbeiPlm}
\docType{class}
\alias{MbeiPlm}


 \title{The MbeiPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class MbeiPlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{MbeiPlm}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{MbeiCnPlm}, \link[aroma.affymetrix]{MbeiSnpPlm}\cr

public static class \bold{MbeiPlm}\cr
extends \link[aroma.affymetrix]{ProbeLevelModel}\cr



  This class represents the multiplicative model in Li \& Wong (2001).
 }

 \usage{MbeiPlm(..., tags="*")}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from ProbeLevelModel}:\cr
calculateResiduals, calculateResidualSet, calculateWeights, fit, getCalculateResidualsFunction, getChipEffects, getChipEffectSet, getProbeAffinities, getProbeAffinityFile, getResidualSet, getWeightsSet

 \bold{Methods inherited from UnitModel}:\cr
findUnitsTodo, getCellIndices, readUnits

 \bold{Methods inherited from Model}:\cr
fit, getAlias, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Model}{
   For a single unit group, the multiplicative model of dChip is:

    \deqn{y_{ij} = \theta_i \phi_j + \varepsilon_{ij}}

   where \eqn{\theta_i} are the chip effects for arrays \eqn{i=1,...,I},
   and \eqn{\phi_j} are the probe affinities for probes \eqn{j=1,...,J}.
   The \eqn{\varepsilon_{ij}} are zero-mean noise with equal variance.
   To make to parameters identifiable, the constraint
   \eqn{\prod_j \phi_j = 1} is added.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
   Li, C. and Wong, W.H. (2001), Genome Biology 2, 1-11.\cr
   Li, C. and Wong, W.H. (2001), Proc. Natl. Acad. Sci USA 98, 31-36.\cr
 }
\keyword{classes}
