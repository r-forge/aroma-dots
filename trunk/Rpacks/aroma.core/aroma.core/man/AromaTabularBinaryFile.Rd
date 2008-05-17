%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaTabularBinaryFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AromaTabularBinaryFile}
\docType{class}
\alias{AromaTabularBinaryFile}


 \title{The AromaTabularBinaryFile class}

 \description{
  Package:  aroma.core \cr
\bold{Class AromaTabularBinaryFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{GenericTabularFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{AromaTabularBinaryFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AromaTabularBinaryFile}\cr
extends \emph{\link[aroma.core]{GenericTabularFile}}\cr



  A AromaTabularBinaryFile represents a file with a binary format.
  It has a well defined header, a data section, and a footer.
 }

 \usage{AromaTabularBinaryFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{GenericTabularFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{[} \tab  -\cr
 \tab \code{[[} \tab  -\cr
 \tab \code{[<-} \tab  -\cr
 \tab \code{\link[aroma.core:allocate.AromaTabularBinaryFile]{allocate}} \tab Creates an AromaTabularBinaryFile.\cr
 \tab \code{as.character} \tab  -\cr
 \tab \code{colMeans} \tab  -\cr
 \tab \code{colMedians} \tab  -\cr
 \tab \code{colStats} \tab  -\cr
 \tab \code{colSums} \tab  -\cr
 \tab \code{dimnames<-} \tab  -\cr
 \tab \code{getBytesPerColumn} \tab  -\cr
 \tab \code{getColClasses} \tab  -\cr
 \tab \code{getColumnNames} \tab  -\cr
 \tab \code{lapply} \tab  -\cr
 \tab \code{nbrOfColumns} \tab  -\cr
 \tab \code{nbrOfRows} \tab  -\cr
 \tab \code{readDataFrame} \tab  -\cr
 \tab \code{\link[aroma.core:readFooter.AromaTabularBinaryFile]{readFooter}} \tab Reads the file footer in XML format into a named nested list.\cr
 \tab \code{subset} \tab  -\cr
 \tab \code{summary} \tab  -\cr
 \tab \code{\link[aroma.core:writeFooter.AromaTabularBinaryFile]{writeFooter}} \tab Writes a named nested list to the file footer in XML format.\cr
}


 \bold{Methods inherited from GenericTabularFile}:\cr
as.character, dim, extractMatrix, getColumnNames, getColumnNameTranslator, nbrOfColumns, nbrOfRows, readColumns, readDataFrame, setColumnNameTranslator, translateColumnNames

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getFullNameTranslator, getName, getPath, getPathname, getTags, gunzip, gzip, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setFullNameTranslator, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link{GenericDataFile}}.
 }
\keyword{classes}
