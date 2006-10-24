%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  UnitGroupsModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{UnitGroupsModel}
\docType{class}
\alias{UnitGroupsModel}


 \title{The UnitGroupsModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class UnitGroupsModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{UnitGroupsModel}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffinePlm}}, \emph{\link[aroma.affymetrix]{AffineSnpPlm}}, \emph{\link[aroma.affymetrix]{MbeiCnPlm}}, \emph{\link[aroma.affymetrix]{MbeiPlm}}, \emph{\link[aroma.affymetrix]{MbeiSnpPlm}}, \emph{\link[aroma.affymetrix]{ProbeLevelModel}}, \emph{\link[aroma.affymetrix]{RmaCnPlm}}, \emph{\link[aroma.affymetrix]{RmaPlm}}, \emph{\link[aroma.affymetrix]{RmaSnpPlm}}\cr

public abstract static class \bold{UnitGroupsModel}\cr
extends \link[R.oo]{Object}\cr



  This class is abstract and represents a generic unit-group model, i.e.
  a model that applies to each group in each unit.  For instance,
  most probeset-summary models such as the RMA model and
  the Li & Wong model belongs to this class of models.
 }

 \usage{UnitGroupsModel(dataSet=NULL, name="modelUnitGroups", path=NULL, ...)}

 \arguments{
   \item{dataSet}{The dataset to which this model should be fitted.}
   \item{path}{The pathname of the directory where to store parameter
     estimates.}
   \item{name}{The name of the unit-group model, which also will be used as
     part of the pathname.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fit.UnitGroupsModel]{fit}} \tab Estimates the model parameters.\cr
 \tab \code{\link[aroma.affymetrix:getCdf.UnitGroupsModel]{getCdf}} \tab Gets the CDF structure for this model.\cr
 \tab \code{\link[aroma.affymetrix:getDataSet.UnitGroupsModel]{getDataSet}} \tab Gets the data set for this model.\cr
 \tab \code{getLabel} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getParameterSet} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getPath.UnitGroupsModel]{getPath}} \tab Gets the path (directory) for this model.\cr
 \tab \code{setLabel} \tab  -\cr
 \tab \code{setup} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson \url{http://www.braju.com/R/}}

 \seealso{
 }
\keyword{classes}
