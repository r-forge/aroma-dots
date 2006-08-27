%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCnRmaModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixCnRmaModel}
\docType{class}
\alias{AffymetrixCnRmaModel}


 \title{The AffymetrixCnRmaModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixCnRmaModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixUnitGroupsModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixRmaModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\emph{\code{AffymetrixCnRmaModel}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{AffymetrixCnRmaModel}\cr
extends \emph{\link[aroma.affymetrix]{AffymetrixRmaModel}}\cr



  This class represents the log-additive model used in RMA.
  It can be used to fit the model on a \code{\link{AffymetrixCelSet}}.
 }

 \usage{AffymetrixCnRmaModel(dataSet=NULL, ..., name="modelCnRma")}

 \arguments{
   \item{dataSet}{An \code{\link{AffymetrixCelSet}} object.}
   \item{path}{The \code{\link[base]{character}} string specifying the path to the directory
      to contain the parameter-estimate files.}
   \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


\bold{Methods inherited from AffymetrixRmaModel}:\cr
getFitFunction, getProbeAffinityClass

\bold{Methods inherited from ProbeLevelModel}:\cr
fit, getChipEffectClass, getChipEffects, getFitFunction, getProbeAffinities, getProbeAffinityClass

\bold{Methods inherited from AffymetrixUnitGroupsModel}:\cr
as.character, fit, getCdf, getDataSet, getParameterSet, getPath, setup

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

   Use \code{\link[aroma.affymetrix:getProbeAffinities.AffymetrixCnRmaModel]{*getProbeAffinities}()} to get the probe-affinity estimates.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \section{Model estimates}{
   The estimated probe affinities are represented by the
   \code{\link{RmaProbeAffinityFile}} class.
 }

 \references{
 }
\keyword{classes}
