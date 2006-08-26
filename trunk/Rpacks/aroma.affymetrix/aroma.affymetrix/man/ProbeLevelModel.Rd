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
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixUnitGroupsModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{ProbeLevelModel}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffymetrixLiWongModel}}, \emph{\link[aroma.affymetrix]{AffymetrixRmaModel}}\cr

public abstract static class \bold{ProbeLevelModel}\cr
extends \emph{\link[aroma.affymetrix]{AffymetrixUnitGroupsModel}}\cr



  This abstract class represents any probe-level model (PLM).
  To quote the \pkg{affyPLM} package:
  "A [...] PLM is a model that is fit to probe-intensity data.
   More specifically, it is where we fit a model with probe level
   and chip level parameters on a probeset by probeset basis",
  where the more general case for a probeset is a unit group
  using Affymetrix CDF terms.
 }

 \usage{ProbeLevelModel(dataSet=NULL, path=filePath("modelPLM", getChipType(getCdf(dataSet))), model=c("pm"), ...)}

 \arguments{
   \item{dataSet}{An \code{\link{AffymetrixCelSet}} object.}
   \item{path}{The \code{\link[base]{character}} string specifying the path to the directory
      to contain the parameter-estimate files.}
   \item{model}{A \code{\link[base]{character}} string specifying how PM and MM values
      should be modelled.  By default only PM signals are used.}
   \item{...}{Arguments passed to \code{\link{AffymetrixUnitGroupsModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:fit.ProbeLevelModel]{fit}} \tab Fits a unit groups (probeset) model.\cr
 \tab \code{\link[aroma.affymetrix:getFitFunction.ProbeLevelModel]{getFitFunction}} \tab Static method to get the low-level function that fits the PLM.\cr
 \tab \code{\link[aroma.affymetrix:getProbeAffinities.ProbeLevelModel]{getProbeAffinities}} \tab Gets the probe affinities for this model.\cr
 \tab \code{\link[aroma.affymetrix:getProbeAffinityClass.ProbeLevelModel]{getProbeAffinityClass}} \tab Static method to get the ProbeAffinity Class object for this model.\cr
}


\bold{Methods inherited from AffymetrixUnitGroupsModel}:\cr
as.character, fit, getCdf, getDataSet, getParameterSet, getPath, setup

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }


 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \section{Model estimates}{
   The estimated probe affinities are represented by the
   \code{\link{ProbeAffinityFile}} class.
   Use \code{\link[aroma.affymetrix:getProbeAffinities.ProbeLevelModel]{*getProbeAffinities}()} to access these.
 }

 \seealso{
   For more details on probe-level models, please see
   the \pkg{affyPLM} package.
 }
\keyword{classes}
