%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaFreqBCghBinaryFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AromaFreqBCghBinaryFile}
\docType{class}
\alias{AromaFreqBCghBinaryFile}


 \title{The AromaFreqBCghBinaryFile class}

 \description{
  Package:  aroma.cn \cr
\bold{Class AromaFreqBCghBinaryFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTabularBinaryFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.cn]{AromaSignalBinaryFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{AromaFreqBCghBinaryFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AromaFreqBCghBinaryFile}\cr
extends \link[aroma.cn]{AromaSignalBinaryFile}\cr



  An AromaFreqBCghBinaryFile is a \code{\link{AromaTabularBinaryFile}}.
 }

 \usage{AromaFreqBCghBinaryFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AromaTabularBinaryFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from AromaSignalBinaryFile}:\cr
allocate, as.character, fromFile, getChipType, getFilenameExtension, getPlatform, nbrOfUnits

 \bold{Methods inherited from AromaTabularBinaryFile}:\cr
[, [[, [<-, allocate, as.character, colMeans, colMedians, colStats, colSums, dim, dimnames<-, extractMatrix, getBytesPerColumn, getColClasses, getColumnNames, lapply, nbrOfColumns, nbrOfRows, readDataFrame, readFooter, readHeader, readRawFooter, setAttributesByTags, subset, summary, updateData, updateDataColumn, writeFooter, writeRawFooter

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
