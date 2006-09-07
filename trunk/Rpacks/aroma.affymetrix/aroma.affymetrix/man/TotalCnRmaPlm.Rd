%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  TotalCnRmaPlm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{TotalCnRmaPlm}
\docType{class}
\alias{TotalCnRmaPlm}


 \title{The TotalCnRmaPlm class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class TotalCnRmaPlm}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{UnitGroupsModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{RmaPlm}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{TotalCnRmaPlm}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{TotalCnRmaPlm}\cr
extends \emph{\link[aroma.affymetrix]{RmaPlm}}\cr



  This class represents the log-additive model used in RMA.
  It can be used to fit the model on a \code{\link{AffymetrixCelSet}}.
 }

 \usage{TotalCnRmaPlm(..., name="modelTotalCnRma")}

 \arguments{
   \item{dataSet}{An \code{\link{AffymetrixCelSet}} object.}
   \item{path}{The \code{\link[base]{character}} string specifying the path to the directory
      to contain the parameter-estimate files.}
   \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


\bold{Methods inherited from RmaPlm}:\cr
getFitFunction, getProbeAffinityClass

\bold{Methods inherited from ProbeLevelModel}:\cr
fit, getChipEffectClass, getChipEffects, getFitFunction, getProbeAffinities, getProbeAffinityClass, highlight, plotMvsPosition

\bold{Methods inherited from UnitGroupsModel}:\cr
as.character, fit, getCdf, getDataSet, getLabel, getName, getParameterSet, getPath, setLabel, setup

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Model}{
   Consider a specific unit group.  The log-additive model in RMA is:

    \deqn{log_2(y_{ij})} = \beta_i + \alpha_j + \eps_{ij}

   where \eqn{\beta_i} are the chip effects for arrays \eqn{i=1,...,I},
   and \eqn{\alpha_j} are the probe affinities for probes \eqn{j=1,...,J}.
   The \eqn{\eps_{ij}} are zero-mean noise with equal variance.

   To minimize the risk for mistakes, all probe-affinity models return
   parameter estimates on the intensity scale.  That is, this class
   returns \eqn{\theta_i = 2^\beta_i} and \eqn{\phi_i = 2^\alpha_i},
   cf. the multiplicative model of Li & Wong.

   Use \code{\link[aroma.affymetrix:getProbeAffinities.TotalCnRmaPlm]{*getProbeAffinities}()} to get the probe-affinity estimates.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \section{Model estimates}{
   The estimated probe affinities are represented by the
   \code{\link{RmaProbeAffinityFile}} class.
 }

 \references{
 }
\keyword{classes}
