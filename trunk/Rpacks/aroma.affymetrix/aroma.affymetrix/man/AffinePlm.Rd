%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffinePlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffinePlm}
\docType{class}
\alias{AffinePlm}


 \title{The AffinePlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffinePlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitGroupsModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\emph{\code{AffinePlm}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffineSnpPlm}}\cr

public abstract static class \bold{AffinePlm}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelModel}}\cr



  This class represents the Bengtsson \& Hössjer (2006) model.
 }

 \usage{AffinePlm(..., background=TRUE, tags="*")}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
   \item{background}{If \code{\link[base:logical]{TRUE}}, background is estimate for each unit group,
     otherwise not. That is, if \code{\link[base:logical]{FALSE}}, a \emph{linear} model without
     offset is fitted, resulting in very similar results as obtained by
     the \code{\link{MbeiPlm}}.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
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
   For a single unit group, the affine model is:

    \deqn{y_{ij} = a + \theta_i \phi_j + \varepsilon_{ij}}

   where \eqn{a} is an offset common to all probe signals,
   \eqn{\theta_i} are the chip effects for arrays \eqn{i=1,...,I},
   and \eqn{\phi_j} are the probe affinities for probes \eqn{j=1,...,J}.
   The \eqn{\varepsilon_{ij}} are zero-mean noise with equal variance.
   The model is constrained such that \eqn{\prod_j \phi_j = 1}.

   Note that with the additional constraint \eqn{a=0} (see arguments above),
   the above model is very similar to \code{\link{MbeiPlm}}.  The differences in
   parameter estimates is due to difference is assumptions about the
   error structure, which in turn affects how the model is estimated.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
   Bengtsson \& Hössjer (2006). \cr
 }
\keyword{classes}
