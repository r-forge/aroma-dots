%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SrmaModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{SrmaModel}
\docType{class}
\alias{SrmaModel}


 \title{The SrmaModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class SrmaModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{ChromosomalModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{SrmaModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{SrmaModel}\cr
extends \emph{\link[aroma.affymetrix]{ChromosomalModel}}\cr



  This class represents the Chromosomal Smoothing Robust Multichip Analysis
  method.
 }

 \usage{SrmaModel(..., bandwidth=10000, tags="*")}

 \arguments{
   \item{cesTuple}{A \code{\link{ChipEffectSetTuple}}.}
   \item{bandwidth}{A single \code{\link[base]{numeric}} specifying the smoothing bandwidth
     in units of nucleotides.}
   \item{...}{Arguments passed to the constructor of
              \code{\link{CopyNumberSegmentationModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{clearCache} \tab  -\cr
 \tab \code{createOutputTuple} \tab  -\cr
 \tab \code{fitOneChromosome} \tab  -\cr
 \tab \code{getBandwidth} \tab  -\cr
 \tab \code{getOutputTuple} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{setBandwidth} \tab  -\cr
}


 \bold{Methods inherited from ChromosomalModel}:\cr
as.character, clearCache, fit, getAlias, getArrays, getAsteriskTag, getChipEffectFiles, getChipType, getChipTypes, getChromosomes, getFullName, getFullNames, getGenome, getGenomeData, getGenomeFile, getListOfChipEffectSets, getListOfGenomeInformations, getName, getNames, getParentPath, getPath, getPcuTheta, getPositionChipTypeUnit, getReportPath, getRootPath, getSetTuple, getTableOfArrays, getTags, nbrOfArrays, nbrOfChipTypes, setAlias, setGenome

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
  \code{\link{CopyNumberSegmentationModel}}.
 }
\keyword{classes}
