%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Model.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{Model}
\docType{class}
\alias{Model}


 \title{The Model class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class Model}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{Model}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffineCnPlm}, \link[aroma.affymetrix]{AffinePlm}, \link[aroma.affymetrix]{AffineSnpPlm}, \link[aroma.affymetrix]{CrlmmModel}, \link[aroma.affymetrix]{MbeiCnPlm}, \link[aroma.affymetrix]{MbeiPlm}, \link[aroma.affymetrix]{MbeiSnpPlm}, \link[aroma.affymetrix]{PlasqModel}, \link[aroma.affymetrix]{ProbeLevelModel}, \link[aroma.affymetrix]{RmaCnPlm}, \link[aroma.affymetrix]{RmaPlm}, \link[aroma.affymetrix]{RmaSnpPlm}, \emph{\link[aroma.affymetrix]{UnitModel}}\cr

public abstract static class \bold{Model}\cr
extends \link[R.oo]{Object}\cr



  This class is abstract and represents a generic model that applies
  to a data set.
 }

 \usage{Model(dataSet=NULL, tags=NULL, ...)}

 \arguments{
   \item{dataSet}{The data set to which this model should be fitted.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be appended to the tags of
      the input data set.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:fit.Model]{fit}} \tab Estimates the model parameters.\cr
 \tab \code{\link[aroma.affymetrix:getDataSet.Model]{getDataSet}} \tab Gets the input data set for this model.\cr
 \tab \code{\link[aroma.affymetrix:getFullName.Model]{getFullName}} \tab Gets the full name of the output set.\cr
 \tab \code{\link[aroma.affymetrix:getName.Model]{getName}} \tab Gets the name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getPath.Model]{getPath}} \tab Gets the path of this model.\cr
 \tab \code{\link[aroma.affymetrix:getRootPath.Model]{getRootPath}} \tab Gets the root path of this model.\cr
 \tab \code{\link[aroma.affymetrix:getTags.Model]{getTags}} \tab Gets the tags of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:setTags.Model]{setTags}} \tab Sets the tags to be appended.\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
 }
\keyword{classes}
