%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FirmaModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{FirmaModel}
\docType{class}
\alias{FirmaModel}


 \title{The FirmaModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class FirmaModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Model}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{FirmaModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{FirmaModel}\cr
extends \emph{\link[aroma.affymetrix]{UnitModel}}\cr



  This class represents the FIRMA (Finding Isoforms using RMA) alternative
  splicing model.

 }

 \usage{FirmaModel(rmaPlm=NULL, summaryMethod="upperQuartile", tags="*", ...)}

 \arguments{
   \item{rmaPlm}{An @RmaPlm object.}
   \item{summaryMethod}{}
   \item{tags}{}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{calculateWeights} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fit.FirmaModel]{fit}} \tab Estimates the model parameters.\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getDataSet} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getFirmaScores.FirmaModel]{getFirmaScores}} \tab Gets the set of FIRMA results for this model.\cr
 \tab \code{\link[aroma.affymetrix:getFitFunction.FirmaModel]{getFitFunction}} \tab Static method to get the low-level function that fits the PLM.\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getPlm} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
}


 \bold{Methods inherited from UnitModel}:\cr
findUnitsTodo, getCellIndices, readUnits

 \bold{Methods inherited from Model}:\cr
fit, getDataSet, getFullName, getName, getPath, getRootPath, getTags, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Ken Simpson (ksimpson[at]wehi.edu.au).}
\keyword{classes}
