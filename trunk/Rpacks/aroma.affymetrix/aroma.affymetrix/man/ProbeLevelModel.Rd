%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeLevelModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ProbeLevelModel}
\docType{class}
\alias{ProbeLevelModel}


 \title{The ProbeLevelModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ProbeLevelModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{UnitGroupsModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{ProbeLevelModel}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffinePlm}}, \emph{\link[aroma.affymetrix]{AffineSnpPlm}}, \emph{\link[aroma.affymetrix]{MbeiCnPlm}}, \emph{\link[aroma.affymetrix]{MbeiPlm}}, \emph{\link[aroma.affymetrix]{MbeiSnpPlm}}, \emph{\link[aroma.affymetrix]{RmaCnPlm}}, \emph{\link[aroma.affymetrix]{RmaPlm}}, \emph{\link[aroma.affymetrix]{RmaSnpPlm}}\cr

public abstract static class \bold{ProbeLevelModel}\cr
extends \emph{\link[aroma.affymetrix]{UnitGroupsModel}}\cr



  This abstract class represents any probe-level model (PLM).
  To quote the \pkg{affyPLM} package:
    "A [...] PLM is a model that is fit to probe-intensity data.
     More specifically, it is where we fit a model with probe level
     and chip level parameters on a probeset by probeset basis",
  where the more general case for a probeset is a unit group
  using Affymetrix CDF terms.
 }

 \usage{ProbeLevelModel(..., name="modelPlm", model=c("pm"), standardize=TRUE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{UnitGroupsModel}}.}
   \item{name}{The name of the model, which is also used in the pathname.}
   \item{model}{A \code{\link[base]{character}} string specifying how PM and MM values
      should be modelled.  By default only PM signals are used.}
   \item{standardize}{If \code{\link[base:logical]{TRUE}}, chip-effect and probe-affinity estimates are
      rescaled such that the product of the probe affinities is one.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:findUnitsTodo.ProbeLevelModel]{findUnitsTodo}} \tab Identifies units for which the PLM has still not been fitted to.\cr
 \tab \code{\link[aroma.affymetrix:fit.ProbeLevelModel]{fit}} \tab Fits a unit groups (probeset) model.\cr
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChipEffects.ProbeLevelModel]{getChipEffects}} \tab Gets the chip-effect files for this model.\cr
 \tab \code{getChipEffectSetClass} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getFitFunction.ProbeLevelModel]{getFitFunction}} \tab Static method to get the low-level function that fits the PLM.\cr
 \tab \code{getFitUnitFunction} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getProbeAffinities.ProbeLevelModel]{getProbeAffinities}} \tab Gets the probe affinities for this model.\cr
 \tab \code{highlight} \tab  -\cr
 \tab \code{plotMvsPosition} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:readUnits.ProbeLevelModel]{readUnits}} \tab Reads data unit by unit.\cr
}


\bold{Methods inherited from UnitGroupsModel}:\cr
as.character, fit, getCdf, getDataSet, getLabel, getName, getParameterSet, getPath, setLabel, setup

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \details{
   In order to minimize the risk for mistakes, but also to be able compare
   results from different PLMs, all PLM subclasses must return parameter
   estimates that meet the following requirements:
   \enumerate{
     \item All parameter estimates must be (stored and returned) on the
       intensity scale, e.g. log-additive models such as \code{\link{RmaPlm}}
       have to transform the parameters on the log-scale to the intensity
       scale.
     \item The probe-affinity estimates \eqn{\phi_j} for a unit group
       must meet the constraint such that \eqn{\prod_j \phi_j = 1},
       or equivalently \eqn{\sum_j \log(\phi_j) = 0}.
   }
   Note that the above probe-affinity constraint guarantees that the
   estimated chip effects across models are on the same scale.
 }

 \author{Henrik Bengtsson \url{http://www.braju.com/R/}}

 \seealso{
   For more details on probe-level models, please see
   the \pkg{affyPLM} package.
 }
\keyword{classes}
