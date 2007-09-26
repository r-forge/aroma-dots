%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CopyNumberSegmentationModel.R
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

 \usage{CopyNumberSegmentationModel(cesTuple=NULL, refTuple=NULL, tags="", genome="Human", ...)}

 \arguments{
   \item{cesTuple}{A \code{\link{ChipEffectSetTuple}}.}
   \item{refTuple}{An optional \code{\link{ChipEffectSetTuple}} for pairwise
     comparisons.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
   \item{genome}{A \code{\link[base]{character}} string specifying what genome is process.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{clearCache} \tab  -\cr
 \tab \code{fit} \tab  -\cr
 \tab \code{fitOne} \tab  -\cr
 \tab \code{getAlias} \tab  -\cr
 \tab \code{getArrays} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getChipTypes} \tab  -\cr
 \tab \code{getFullNames} \tab  -\cr
 \tab \code{getListOfChipEffectSets} \tab  -\cr
 \tab \code{getListOfReferences} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getNames} \tab  -\cr
 \tab \code{getReferenceFiles} \tab  -\cr
 \tab \code{getRegions} \tab  -\cr
 \tab \code{getTableOfArrays} \tab  -\cr
 \tab \code{migrateTool} \tab  -\cr
 \tab \code{nbrOfArrays} \tab  -\cr
 \tab \code{nbrOfChipTypes} \tab  -\cr
 \tab \code{plot} \tab  -\cr
 \tab \code{setAlias} \tab  -\cr
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
