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
\code{~~+--}\code{\link[aroma.affymetrix]{UnitGroupsModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{MbeiPlm}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{MbeiCnPlm}}, \emph{\link[aroma.affymetrix]{MbeiSnpPlm}}\cr

public abstract static class \bold{MbeiPlm}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelModel}}\cr



  This class represents the Li \& Wong (2001) model.
 }

 \usage{MbeiPlm(..., name="modelMbeiPlm")}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
   \item{name}{The name of the model, which is also used in the pathname.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getProbeAffinities} \tab  -\cr
}


\bold{Methods inherited from ProbeLevelModel}:\cr
findUnitsTodo, fit, getCellIndices, getChipEffects, getChipEffectSetClass, getFitFunction, getFitUnitFunction, getProbeAffinities, highlight, plotMvsPosition, readUnits

\bold{Methods inherited from UnitGroupsModel}:\cr
as.character, fit, getCdf, getDataSet, getLabel, getName, getParameterSet, getPath, setLabel, setup

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Model}{
   For a single unit group, the multiplicative model of dChip is:

    \deqn{y_{ij} = \theta_i \phi_j + \varepsilon_{ij}}

   where \eqn{\theta_i} are the chip effects for arrays \eqn{i=1,...,I},
   and \eqn{\phi_j} are the probe affinities for probes \eqn{j=1,...,J}.
   The \eqn{\varepsilon_{ij}} are zero-mean noise with equal variance.

   In addition, we modify th constraint such that it is guaranteed that
   \eqn{\prod_j \phi_j = 1}.
 }

 \author{Henrik Bengtsson \url{http://www.braju.com/R/}}

 \references{
   Li, C. and Wong, W.H. (2001), Genome Biology 2, 1-11.\cr
   Li, C. and Wong, W.H. (2001), Proc. Natl. Acad. Sci USA 98, 31-36.\cr
 }
\keyword{classes}
