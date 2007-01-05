%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PlasqModel.R
% 
% on Sat Jan 06 01:54:45 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{PlasqModel}
\docType{class}
\alias{PlasqModel}


 \title{The PlasqModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class PlasqModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{UnitModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{PlasqModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{PlasqModel}\cr
extends \emph{\link[aroma.affymetrix]{UnitModel}}\cr



  This class represents the PLASQ model.
 }

 \usage{PlasqModel(...)}

 \arguments{
   \item{...}{Arguments passed to the constructor of \code{\link{UnitModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:findUnitsTodo.PlasqModel]{findUnitsTodo}} \tab Identifies units for which the PLM has still not been fitted to.\cr
 \tab \code{\link[aroma.affymetrix:fit.PlasqModel]{fit}} \tab Fits the model.\cr
 \tab \code{getChipEffects} \tab  -\cr
 \tab \code{getProbeParameters} \tab  -\cr
}


\bold{Methods inherited from UnitModel}:\cr
as.character, findUnitsTodo, fit, getCdf, getCellIndices, getDataSet, getFitUnitFunction, getFullName, getName, getParameterSet, getPath, getRootPath, getTags, readUnits, setTags

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getInternalAddress, getStaticInstance, hasField, hashCode, isReferable, ll, load, names, newInstance, novirtual, objectSize, print, save, staticCode


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
