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
\code{~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{UnitGroupsModel}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffinePlm}}, \emph{\link[aroma.affymetrix]{AffineSnpPlm}}, \emph{\link[aroma.affymetrix]{MbeiCnPlm}}, \emph{\link[aroma.affymetrix]{MbeiPlm}}, \emph{\link[aroma.affymetrix]{MbeiSnpPlm}}, \emph{\link[aroma.affymetrix]{ProbeLevelModel}}, \emph{\link[aroma.affymetrix]{RmaCnPlm}}, \emph{\link[aroma.affymetrix]{RmaPlm}}, \emph{\link[aroma.affymetrix]{RmaSnpPlm}}\cr

public abstract static class \bold{UnitGroupsModel}\cr
extends \emph{\link[aroma.affymetrix]{UnitModel}}\cr



  This class is abstract and represents a generic unit-group model, i.e.
  a model that applies to each group in each unit.  For instance,
  most probeset-summary models such as the RMA model and
  the Li & Wong model belongs to this class of models.
 }

 \usage{UnitGroupsModel(...)}

 \arguments{
   \item{...}{Arguments passed to constructor of superclass \code{\link{UnitModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getLabel} \tab  -\cr
 \tab \code{setLabel} \tab  -\cr
 \tab \code{setup} \tab  -\cr
}


\bold{Methods inherited from UnitModel}:\cr
as.character, findUnitsTodo, fit, getCdf, getCellIndices, getDataSet, getFitUnitFunction, getFullName, getName, getParameterSet, getPath, getRootPath, getTags, readUnits, setTags

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
 }
\keyword{classes}
