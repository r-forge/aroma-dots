%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipGenomeInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{DChipGenomeInformation}
\docType{class}
\alias{DChipGenomeInformation}


 \title{The DChipGenomeInformation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class DChipGenomeInformation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{GenomeInformation}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{DChipGenomeInformation}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{DChipGenomeInformation}\cr
extends \emph{\link[aroma.affymetrix]{GenomeInformation}}\cr



  This class represents dChip genome information files, which typically
  contains information about chromosomal locations of the units.
 }

 \usage{DChipGenomeInformation(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{GenomeInformation}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:fromChipType.DChipGenomeInformation]{fromChipType}} \tab Defines a DChipGenomeInformation object by chip type.\cr
 \tab \code{readData} \tab  -\cr
}


 \bold{Methods inherited from GenomeInformation}:\cr
fromCdf, fromChipType, fromDataSet, getChipType, getChromosomes, getChromosomeStats, getData, getPositions, getUnitIndices, getUnitsOnChromosome, nbrOfUnits, plotDensity, readData

 \bold{Methods inherited from AffymetrixFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, readChecksum, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \details{
   The dChip genome information files for various chip types can be
   downloaded from \url{http://www.dchip.org/}.  Put each file in a
   directory named identically as the corresponding chip type under the
   \emph{annotations/} directory, e.g.
   \emph{annotations/Mapping50K\_Hind240/50k hind genome info AfAm
   june 05 hg17.xls}.
   Note that dChip changes the filename and file format slightly between
   chip types, but currently the \code{\link[aroma.affymetrix:fromChipType.DChipGenomeInformation]{*fromChipType}()} basically searches
   for files with names consisting of \code{"genome info"} or
   \code{"genome_info"}.  At least for the most common chip types, there
   is no need to rename the files in order for this class to recognize them.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
