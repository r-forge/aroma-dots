%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaUnitTabularBinaryFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AromaUnitTabularBinaryFile}
\docType{class}
\alias{AromaUnitTabularBinaryFile}


 \title{The AromaUnitTabularBinaryFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AromaUnitTabularBinaryFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTabularBinaryFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{AromaUnitTabularBinaryFile}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AromaUfcFile}}, \emph{\link[aroma.affymetrix]{AromaUflFile}}, \link[aroma.affymetrix]{AromaUgpFile}\cr

public abstract static class \bold{AromaUnitTabularBinaryFile}\cr
extends \link[aroma.core]{AromaTabularBinaryFile}\cr



  A AromaUnitTabularBinaryFile is an \code{\link{AromaTabularBinaryFile}} with
  the constraint that the rows map one-to-one to, and in the same order as,
  the units in a CDF.  The (full) chip type of the CDF is given by the
  mandatory file footer \code{chipType}.
 }

 \usage{AromaUnitTabularBinaryFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AromaTabularBinaryFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{allocateFromCdf} \tab  -\cr
 \tab \code{as.character} \tab  -\cr
 \tab \code{byChipType} \tab  -\cr
 \tab \code{fromChipType} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getFilenameExtension} \tab  -\cr
 \tab \code{getImage} \tab  -\cr
 \tab \code{importFrom} \tab  -\cr
 \tab \code{importFromAffymetrixTabularFile} \tab  -\cr
 \tab \code{importFromDChipGenomeInformation} \tab  -\cr
 \tab \code{importFromGenericTabularFile} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
}


 \bold{Methods inherited from AromaTabularBinaryFile}:\cr
[, [[, [<-, allocate, as.character, colMeans, colMedians, colStats, colSums, dim, dimnames<-, getBytesPerColumn, getColClasses, getColumnNames, lapply, nbrOfColumns, nbrOfRows, readDataFrame, readFooter, readHeader, readRawFooter, setAttributesByTags, subset, summary, updateData, updateDataColumn, writeFooter, writeRawFooter

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link{AromaTabularBinaryFile}}.
 }
\keyword{classes}
