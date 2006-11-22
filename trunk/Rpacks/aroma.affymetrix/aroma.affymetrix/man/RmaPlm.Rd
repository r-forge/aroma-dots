%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RmaPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{RmaPlm}
\docType{class}
\alias{RmaPlm}


 \title{The RmaPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class RmaPlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{UnitGroupsModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{RmaPlm}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{RmaCnPlm}}, \emph{\link[aroma.affymetrix]{RmaSnpPlm}}\cr

public abstract static class \bold{RmaPlm}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelModel}}\cr



  This class represents the log-additive model used in RMA.
 }

 \usage{RmaPlm(..., tags="*")}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getProbeAffinities} \tab  -\cr
}


\bold{Methods inherited from ProbeLevelModel}:\cr
findUnitsTodo, fit, getCellIndices, getChipEffects, getChipEffectSetClass, getFitFunction, getFitUnitFunction, getProbeAffinities, highlight, plotMvsPosition, readUnits

\bold{Methods inherited from UnitGroupsModel}:\cr
as.character, fit, getCdf, getDataSet, getFullName, getLabel, getName, getParameterSet, getPath, getRootPath, getTags, setLabel, setTags, setup

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Model}{
   For a single unit group, the log-additive model of RMA is:

    \deqn{log_2(y_{ij}) = \beta_i + \alpha_j + \varepsilon_{ij}}

   where \eqn{\beta_i} are the chip effects for arrays \eqn{i=1,...,I},
   and \eqn{\alpha_j} are the probe affinities for probes \eqn{j=1,...,J}.
   The \eqn{\varepsilon_{ij}} are zero-mean noise with equal variance.
   The model is constrained such that \eqn{\sum_j{\alpha_j} = 0}.

   Note that all PLM classes must return parameters on the intensity scale.
   For this class that means that \eqn{\theta_i = 2^\beta_i} and
   \eqn{\phi_i = 2^\alpha_i} are returned.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
  Irizarry et al. \emph{Summaries of Affymetrix GeneChip probe level data}.
  NAR, 2003, 31, e15.\cr
 }
\keyword{classes}
