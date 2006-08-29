%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixTotalCnLiWongModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixTotalCnLiWongModel}
\docType{class}
\alias{AffymetrixTotalCnLiWongModel}


 \title{The AffymetrixTotalCnLiWongModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixTotalCnLiWongModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixUnitGroupsModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixLiWongModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{AffymetrixTotalCnLiWongModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AffymetrixTotalCnLiWongModel}\cr
extends \emph{\link[aroma.affymetrix]{AffymetrixLiWongModel}}\cr



  This class represents the multiplicate model by Li and Wong (2001).
  It can be used to fit the model on a \code{\link{AffymetrixCelSet}}.
 }

 \usage{AffymetrixTotalCnLiWongModel(..., name="modelTotalCnLiWong")}

 \arguments{
   \item{dataSet}{An \code{\link{AffymetrixCelSet}} object.}
   \item{path}{The \code{\link[base]{character}} string specifying the path to the directory
      to contain the parameter-estimate files.}
   \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


\bold{Methods inherited from AffymetrixLiWongModel}:\cr
getChipEffects2, getFitFunction, getProbeAffinityClass

\bold{Methods inherited from ProbeLevelModel}:\cr
fit, getChipEffectClass, getChipEffects, getFitFunction, getProbeAffinities, getProbeAffinityClass, highlight, plotMvsPosition

\bold{Methods inherited from AffymetrixUnitGroupsModel}:\cr
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

   Use \code{\link[aroma.affymetrix:getProbeAffinities.AffymetrixTotalCnLiWongModel]{*getProbeAffinities}()} to get the probe-affinity estimates.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \section{Model estimates}{
   The estimated probe affinities are represented by the
   \code{\link{LiWongProbeAffinityFile}} class.
 }

 \references{
 }
\keyword{classes}
