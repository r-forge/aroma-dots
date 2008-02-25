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
  Package:  aroma.affymetrix \cr
\bold{Class AromaTabularBinaryFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{AromaTabularBinaryFile}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AromaUfcFile}}, \emph{\link[aroma.affymetrix]{AromaUflFile}}, \link[aroma.affymetrix]{AromaUgpFile}, \emph{\link[aroma.affymetrix]{AromaUnitTabularBinaryFile}}\cr

public static class \bold{AromaTabularBinaryFile}\cr
extends \link[aroma.affymetrix]{GenericDataFile}\cr



  A AromaTabularBinaryFile represents a file with a binary format.
  It has a well defined header, a data section, and a footer.
 }

 \usage{AromaTabularBinaryFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{GenericDataFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{[} \tab  -\cr
 \tab \code{[[} \tab  -\cr
 \tab \code{[<-} \tab  -\cr
 \tab \code{allocate} \tab  -\cr
 \tab \code{as.character} \tab  -\cr
 \tab \code{colMeans} \tab  -\cr
 \tab \code{colMedians} \tab  -\cr
 \tab \code{colnames} \tab  -\cr
 \tab \code{colStats} \tab  -\cr
 \tab \code{colSums} \tab  -\cr
 \tab \code{dim} \tab  -\cr
 \tab \code{dimnames<-} \tab  -\cr
 \tab \code{getBytesPerColumn} \tab  -\cr
 \tab \code{getColClasses} \tab  -\cr
 \tab \code{lapply} \tab  -\cr
 \tab \code{nbrOfColumns} \tab  -\cr
 \tab \code{nbrOfRows} \tab  -\cr
 \tab \code{ncol} \tab  -\cr
 \tab \code{nrow} \tab  -\cr
 \tab \code{readFooter} \tab  -\cr
 \tab \code{subset} \tab  -\cr
 \tab \code{summary} \tab  -\cr
 \tab \code{writeFooter} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link{GenericDataFile}}.
 }
\keyword{classes}
