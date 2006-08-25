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
\link[aroma.affymetrix]{AffymetrixLiWongModel}\cr

public abstract static class \bold{AffymetrixUnitGroupsModel}\cr
extends \link[R.oo]{Object}\cr



  This class is abstract and represents a generic unit-group model, i.e.
  a model that applies to each group in each unit.  For instance,
  most probeset-summary models such as the RMA model and
  the Li & Wong model belongs to this class of models.
 }

 \usage{AffymetrixUnitGroupsModel(path=NULL, ..., dataSet=NULL)}

 \arguments{
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{fit} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getDataSet} \tab  -\cr
 \tab \code{getParameterSet} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{setup} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
 }
\keyword{classes}
