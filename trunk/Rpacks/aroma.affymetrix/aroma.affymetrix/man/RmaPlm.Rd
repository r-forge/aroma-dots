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
\code{~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitGroupsModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\emph{\code{RmaPlm}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{RmaCnPlm}}, \emph{\link[aroma.affymetrix]{RmaSnpPlm}}\cr

public abstract static class \bold{RmaPlm}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelModel}}\cr



  This class represents the log-additive model used in RMA.
 }

 \usage{RmaPlm(..., tags="*", flavor=c("affyPLM", "affyPLMold", "oligo"))}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
   \item{flavor}{A \code{\link[base]{character}} string specifying what model fitting algorithm
     to be used.  This makes it possible to get identical estimates as other
     packages.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getProbeAffinities} \tab  -\cr
}


\bold{Methods inherited from ProbeLevelModel}:\cr
findUnitsTodo, fit, getCellIndices, getChipEffects, getChipEffectSetClass, getFitFunction, getFitUnitFunction, getProbeAffinities, getRootPath, highlight, plotMvsPosition, readUnits

\bold{Methods inherited from UnitGroupsModel}:\cr
getLabel, setLabel, setup

\bold{Methods inherited from UnitModel}:\cr
as.character, findUnitsTodo, fit, getCdf, getCellIndices, getDataSet, getFitUnitFunction, getFullName, getName, getParameterSet, getPath, getRootPath, getTags, readUnits, setTags

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

 \section{Different flavors of model fitting}{
   There are a few differ algorithms available for fitting the same
   probe-level model.  The default method (\code{flavor="affyPLM"})
   implements in the \pkg{affyPLM} package fits the model parameters
   robustly using an M-estimator.

   In addition to the above method, which is the preferred one, additional
   algorithms have been ported/interfaced to.  Note that these are often
   provided for the purpose of replicating the results of other packages,
   but they might not be complete in the sense of what a
   \code{\link{ProbeLevelModel}} should return or how the parameters are
   constrained.
   The algorithm (\code{flavor="oligo"}) used by the \pkg{oligo} package,
   which originates from the \pkg{affy} packages, fits the model using
   median polish, which is a non-robust estimator.  Note that this algorithm
   does not constraint the probe-effect parameters to multiply to one on
   the intensity scale. As a matter of fact, the probe-effect parameters
   are not even returned (defaulting to ones).
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
  Irizarry et al. \emph{Summaries of Affymetrix GeneChip probe level data}.
  NAR, 2003, 31, e15.\cr
 }
\keyword{classes}
