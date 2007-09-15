%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CopyNumberSegmentationModel.R
% 
% on Sat Sep 15 13:21:17 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CopyNumberSegmentationModel}
\docType{class}
\alias{CopyNumberSegmentationModel}


 \title{The CopyNumberSegmentationModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CopyNumberSegmentationModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{CopyNumberSegmentationModel}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{CbsModel}, \link[aroma.affymetrix]{GladModel}\cr

public abstract static class \bold{CopyNumberSegmentationModel}\cr
extends \link[R.oo]{Object}\cr



  This \emph{abstract} class represents a copy-number segmentation model.
 }

 \usage{CopyNumberSegmentationModel(cesTuple=NULL, referenceList=NULL, tags="", genome="Human", ...)}

 \arguments{
   \item{cesTuple}{A \code{\link{ChipEffectSetTuple}}.}
   \item{referenceList}{A single or a \code{\link[base]{list}} of \code{\link{ChipEffectFile}}:s.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
   \item{genome}{A \code{\link[base]{character}} string specifying what genome is process.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:fit.CopyNumberSegmentationModel]{fit}} \tab Fits the model.\cr
 \tab \code{\link[aroma.affymetrix:fitOne.CopyNumberSegmentationModel]{fitOne}} \tab Fits the copy-number segmentation model for one chromosome in one sample.\cr
 \tab \code{getAlias} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getArrays.CopyNumberSegmentationModel]{getArrays}} \tab Gets the names of the arrays.\cr
 \tab \code{getChipEffectFiles} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChipType.CopyNumberSegmentationModel]{getChipType}} \tab Gets a label for all chip types merged.\cr
 \tab \code{getChipTypes} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChromosomes.CopyNumberSegmentationModel]{getChromosomes}} \tab Gets the chromosomes available.\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getFullNames} \tab  -\cr
 \tab \code{getGenome} \tab  -\cr
 \tab \code{getListOfChipEffects} \tab  -\cr
 \tab \code{getListOfGenomeInformations} \tab  -\cr
 \tab \code{getListOfReferences} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getNames} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getReferenceFiles} \tab  -\cr
 \tab \code{getRegions} \tab  -\cr
 \tab \code{getSetTuple} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getTableOfArrays.CopyNumberSegmentationModel]{getTableOfArrays}} \tab Gets a table of arrays.\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{migrateTool} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:nbrOfArrays.CopyNumberSegmentationModel]{nbrOfArrays}} \tab Gets the number of arrays.\cr
 \tab \code{\link[aroma.affymetrix:nbrOfChipTypes.CopyNumberSegmentationModel]{nbrOfChipTypes}} \tab Gets the number of chip types.\cr
 \tab \code{plot} \tab  -\cr
 \tab \code{setAlias} \tab  -\cr
 \tab \code{setGenome} \tab  -\cr
 \tab \code{writeRegions} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Requirements}{
   This class requires genome information annotation files for
   every chip type.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
