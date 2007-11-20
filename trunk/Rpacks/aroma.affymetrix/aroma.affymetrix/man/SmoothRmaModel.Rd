%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SmoothSaModel.R
% 
% on Tue Nov 20 12:40:57 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{SmoothRmaModel}
\docType{class}
\alias{SmoothRmaModel}


 \title{The SmoothRmaModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class SmoothRmaModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{ChromosomalModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{SmoothMultiarrayModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{SmoothRmaModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{SmoothRmaModel}\cr
extends \emph{\link[aroma.affymetrix]{SmoothMultiarrayModel}}\cr



  This class represents the Chromosomal Smoothing Robust Multichip Analysis
  method.
 }

 \usage{SmoothRmaModel(...)}

 \arguments{
   \item{...}{Arguments passed to the constructor of
              \code{\link{SmoothMultiarrayModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getFitFunction} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
}


 \bold{Methods inherited from SmoothMultiarrayModel}:\cr
as.character, clearCache, createOutputTuple, fitOneChromosome, getAsteriskTag, getBandwidth, getFitFunction, getOutputTuple, getRootPath, setBandwidth

 \bold{Methods inherited from ChromosomalModel}:\cr
as.character, clearCache, fit, getAlias, getArrays, getAsteriskTag, getChipEffectFiles, getChipType, getChipTypes, getChromosomes, getFullName, getFullNames, getGenome, getGenomeData, getGenomeFile, getListOfChipEffectSets, getListOfGenomeInformations, getName, getNames, getParentPath, getPath, getPcuTheta, getPositionChipTypeUnit, getReportPath, getRootPath, getSetTuple, getTableOfArrays, getTags, nbrOfArrays, nbrOfChipTypes, setAlias, setGenome

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
