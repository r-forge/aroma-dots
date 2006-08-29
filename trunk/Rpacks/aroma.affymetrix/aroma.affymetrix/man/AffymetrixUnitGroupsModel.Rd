%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixUnitGroupsModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixUnitGroupsModel}
\docType{class}
\alias{AffymetrixUnitGroupsModel}


 \title{The AffymetrixUnitGroupsModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixUnitGroupsModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{AffymetrixUnitGroupsModel}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffymetrixAffineModel}}, \emph{\link[aroma.affymetrix]{AffymetrixLiWongModel}}, \emph{\link[aroma.affymetrix]{AffymetrixRmaModel}}, \link[aroma.affymetrix]{AffymetrixTotalCnAffineModel}, \link[aroma.affymetrix]{AffymetrixTotalCnLiWongModel}, \link[aroma.affymetrix]{AffymetrixTotalCnRmaModel}, \emph{\link[aroma.affymetrix]{ProbeLevelModel}}\cr

public abstract static class \bold{AffymetrixUnitGroupsModel}\cr
extends \link[R.oo]{Object}\cr



  This class is abstract and represents a generic unit-group model, i.e.
  a model that applies to each group in each unit.  For instance,
  most probeset-summary models such as the RMA model and
  the Li & Wong model belongs to this class of models.
 }

 \usage{AffymetrixUnitGroupsModel(dataSet=NULL, name="modelUnitGroups", path=NULL, ...)}

 \arguments{
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fit.AffymetrixUnitGroupsModel]{fit}} \tab Estimates the model parameters.\cr
 \tab \code{\link[aroma.affymetrix:getCdf.AffymetrixUnitGroupsModel]{getCdf}} \tab Gets the CDF structure for this model.\cr
 \tab \code{\link[aroma.affymetrix:getDataSet.AffymetrixUnitGroupsModel]{getDataSet}} \tab Gets the data set for this model.\cr
 \tab \code{getLabel} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getParameterSet} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getPath.AffymetrixUnitGroupsModel]{getPath}} \tab Gets the path (directory) for this model.\cr
 \tab \code{setLabel} \tab  -\cr
 \tab \code{setup} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
 }
\keyword{classes}
