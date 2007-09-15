%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CbsModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CbsModel}
\docType{class}
\alias{CbsModel}


 \title{The CbsModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CbsModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{CopyNumberSegmentationModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{CbsModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CbsModel}\cr
extends \emph{\link[aroma.affymetrix]{CopyNumberSegmentationModel}}\cr



  This class represents the Circular Binary Segmentation (CBS) model [1].
 }

 \usage{CbsModel(cesTuple=NULL, ...)}

 \arguments{
   \item{cesTuple}{A \code{\link{ChipEffectSetTuple}}.}
   \item{...}{Arguments passed to the constructor of
              \code{\link{CopyNumberSegmentationModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from CopyNumberSegmentationModel}:\cr
as.character, calculateChromosomeStatistics, fit, fitOne, getAlias, getArrays, getAsteriskTag, getChipEffectFiles, getChipType, getChipTypes, getChromosomes, getFullName, getFullNames, getGenome, getGenomeData, getGenomeFile, getListOfChipEffects, getListOfGenomeInformations, getListOfReferences, getName, getNames, getPath, getRawCnData, getReferenceFiles, getRegions, getReportPath, getSetTuple, getTableOfArrays, getTags, migrateTool, nbrOfArrays, nbrOfChipTypes, newPlot, plot, setAlias, setGenome, writeRegions

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
  \code{\link{CopyNumberSegmentationModel}}.
 }

 \references{
  [1] Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M.
      \emph{Circular binary segmentation for the analysis of array-based
      DNA copy number data. Biostatistics 5: 557-572, 2004.}\cr
  [2] Venkatraman, E. S. & Olshen, A. B.
      \emph{A faster circular binary segmentation algorithm for the
      analysis of array CGH data}. Bioinformatics, 2007.\cr
 }
\keyword{classes}
