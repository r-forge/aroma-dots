%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipDcpFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{DChipDcpFile}
\docType{class}
\alias{DChipDcpFile}


 \title{The DChipDcpFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class DChipDcpFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{DChipDcpFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{DChipDcpFile}\cr
extends \link[aroma.affymetrix]{AffymetrixFile}\cr



  A DChipDcpFile object represents a DChip DCP file.
 }

 \usage{DChipDcpFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{dim} \tab  -\cr
 \tab \code{fromFile} \tab  -\cr
 \tab \code{getCalls} \tab  -\cr
 \tab \code{getExcludes} \tab  -\cr
 \tab \code{getFileFormat} \tab  -\cr
 \tab \code{getHeader} \tab  -\cr
 \tab \code{getNormalizedIntensities} \tab  -\cr
 \tab \code{getRawIntensities} \tab  -\cr
 \tab \code{getThetas} \tab  -\cr
 \tab \code{getThetasAB} \tab  -\cr
 \tab \code{getThetaStds} \tab  -\cr
 \tab \code{hasMbeiData} \tab  -\cr
 \tab \code{hasNormalizedData} \tab  -\cr
 \tab \code{nbrOfCells} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
}


 \bold{Methods inherited from AffymetrixFile}:\cr
getAttribute, getAttributes, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, testAttributes

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \seealso{
   \code{\link{DChipDcpSet}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
