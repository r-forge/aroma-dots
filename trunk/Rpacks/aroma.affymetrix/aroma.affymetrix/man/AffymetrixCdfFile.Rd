%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixCdfFile}
\docType{class}
\alias{AffymetrixCdfFile}


 \title{The AffymetrixCdfFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixCdfFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaMicroarrayDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AromaChipTypeAnnotationFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{UnitNamesFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{AffymetrixCdfFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AffymetrixCdfFile}\cr
extends \emph{\link[aroma.core]{UnitNamesFile}}\cr



  An AffymetrixCdfFile object represents a generic Affymetrix CDF file.
 }

 \usage{AffymetrixCdfFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AromaChipTypeAnnotationFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{convert} \tab  -\cr
 \tab \code{createMonoCell} \tab  -\cr
 \tab \code{getCellIndices} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getDimension} \tab  -\cr
 \tab \code{getFileFormat} \tab  -\cr
 \tab \code{getGenomeInformation} \tab  -\cr
 \tab \code{getImage} \tab  -\cr
 \tab \code{getMonoCell} \tab  -\cr
 \tab \code{getPlatform} \tab  -\cr
 \tab \code{getUnitGroupNamesFromUgcMap} \tab  -\cr
 \tab \code{getUnitNames} \tab  -\cr
 \tab \code{getUnitTypes} \tab  -\cr
 \tab \code{hasUnitTypes} \tab  -\cr
 \tab \code{isMonocellCdf} \tab  -\cr
 \tab \code{isPm} \tab  -\cr
 \tab \code{nbrOfCells} \tab  -\cr
 \tab \code{nbrOfColumns} \tab  -\cr
 \tab \code{nbrOfQcUnits} \tab  -\cr
 \tab \code{nbrOfRows} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{readDataFrame} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
 \tab \code{setGenomeInformation} \tab  -\cr
 \tab \code{setSnpInformation} \tab  -\cr
}


 \bold{Methods inherited from UnitNamesFile}:\cr
getChipType, getPlatform, getUnitNames, indexOf, nbrOfUnits

 \bold{Methods inherited from AromaChipTypeAnnotationFile}:\cr
byChipType, byName, findByChipType, fromChipType, fromFile, getChipType, getHeader, getPlatform

 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getPlatform

 \bold{Methods inherited from AromaMicroarrayDataFile}:\cr
getAttributeXY, getChipType, getPlatform, getPloidy, hasAttributeXY, setAttributesByTags, setAttributeXY

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getFullNameTranslator, getName, getPath, getPathname, getTags, gunzip, gzip, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setFullName, setFullNameTranslator, setName, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
