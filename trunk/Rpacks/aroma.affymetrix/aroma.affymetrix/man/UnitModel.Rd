%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  UnitModel.R
% 
% on Sat Jan 06 01:58:47 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{UnitModel}
\docType{class}
\alias{UnitModel}


 \title{The UnitModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class UnitModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{UnitModel}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffinePlm}}, \emph{\link[aroma.affymetrix]{AffineSnpPlm}}, \emph{\link[aroma.affymetrix]{MbeiCnPlm}}, \emph{\link[aroma.affymetrix]{MbeiPlm}}, \emph{\link[aroma.affymetrix]{MbeiSnpPlm}}, \link[aroma.affymetrix]{PlasqModel}, \emph{\link[aroma.affymetrix]{ProbeLevelModel}}, \emph{\link[aroma.affymetrix]{RmaCnPlm}}, \emph{\link[aroma.affymetrix]{RmaPlm}}, \emph{\link[aroma.affymetrix]{RmaSnpPlm}}, \emph{\link[aroma.affymetrix]{UnitGroupsModel}}\cr

public abstract static class \bold{UnitModel}\cr
extends \link[R.oo]{Object}\cr



  This class is abstract and represents a generic unit model, i.e.
  a model that is applied to each unit separately.
 }

 \usage{UnitModel(dataSet=NULL, tags="", ...)}

 \arguments{
   \item{dataSet}{The data set to which this model should be fitted.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be appended to the tags of
      the input data set.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{findUnitsTodo} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fit.UnitModel]{fit}} \tab Estimates the model parameters.\cr
 \tab \code{\link[aroma.affymetrix:getCdf.UnitModel]{getCdf}} \tab Gets the CDF structure for this model.\cr
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getDataSet.UnitModel]{getDataSet}} \tab Gets the data set for this model.\cr
 \tab \code{getFitUnitFunction} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getFullName.UnitModel]{getFullName}} \tab Gets the full name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getName.UnitModel]{getName}} \tab Gets the name of the output summarized data set.\cr
 \tab \code{getParameterSet} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getPath.UnitModel]{getPath}} \tab Gets the path of this model.\cr
 \tab \code{\link[aroma.affymetrix:getRootPath.UnitModel]{getRootPath}} \tab Gets the root path of this model.\cr
 \tab \code{\link[aroma.affymetrix:getTags.UnitModel]{getTags}} \tab Gets the tags of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:readUnits.UnitModel]{readUnits}} \tab Reads data unit by unit.\cr
 \tab \code{setTags} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getInternalAddress, getStaticInstance, hasField, hashCode, isReferable, ll, load, names, newInstance, novirtual, objectSize, print, save, staticCode


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
 }
\keyword{classes}
