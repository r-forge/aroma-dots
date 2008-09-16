%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaTotalCnBinaryFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AromaTotalCnBinaryFile}
\docType{class}
\alias{AromaTotalCnBinaryFile}


 \title{The AromaTotalCnBinaryFile class}

 \description{
  Package:  aroma.cn \cr
\bold{Class AromaTotalCnBinaryFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{GenericTabularFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaTabularBinaryFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaPlatformInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.cn]{AromaSignalBinaryFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~~~~~~+--}\code{AromaTotalCnBinaryFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AromaTotalCnBinaryFile}\cr
extends \link[aroma.cn]{AromaSignalBinaryFile}\cr



  An AromaTotalCnBinaryFile is a \code{\link[aroma.core]{AromaTabularBinaryFile}}.
 }

 \usage{AromaTotalCnBinaryFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link[aroma.core]{AromaTabularBinaryFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{extractRawCopyNumbers} \tab  -\cr
}


 \bold{Methods inherited from AromaSignalBinaryFile}:\cr
allocate, as.character, fromFile, getAromaUgpFile, getChipType, getFilenameExtension, getPlatform, getUnitNamesFile, nbrOfUnits

 \bold{Methods inherited from AromaPlatformInterface}:\cr
getAromaPlatform, getPlatform

 \bold{Methods inherited from AromaTabularBinaryFile}:\cr
[, [[, [<-, allocate, as.character, colMeans, colMedians, colStats, colSums, dimnames<-, getBytesPerColumn, getColClasses, getColumnNames, importFrom, lapply, nbrOfColumns, nbrOfRows, readColumns, readDataFrame, readFooter, readHeader, readRawFooter, setAttributesByTags, subset, summary, updateData, updateDataColumn, writeFooter, writeRawFooter

 \bold{Methods inherited from GenericTabularFile}:\cr
as.character, dim, extractMatrix, getColumnNames, getColumnNameTranslator, nbrOfColumns, nbrOfRows, readColumns, readDataFrame, setColumnNameTranslator, translateColumnNames

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getCreatedOn, getFilename, getFileSize, getFileType, getFullName, getFullNameTranslator, getLastAccessedOn, getLastModifiedOn, getName, getPath, getPathname, getTags, gunzip, gzip, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setFullName, setFullNameTranslator, setName, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
