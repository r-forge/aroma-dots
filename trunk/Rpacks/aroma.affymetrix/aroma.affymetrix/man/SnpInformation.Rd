%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{SnpInformation}
\docType{class}
\alias{SnpInformation}


 \title{The SnpInformation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class SnpInformation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{SnpInformation}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{DChipSnpInformation}, \link[aroma.affymetrix]{UflSnpInformation}\cr

public abstract static class \bold{SnpInformation}\cr
extends \link[aroma.affymetrix]{GenericDataFile}\cr


 }

 \usage{SnpInformation(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{GenericDataFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fromCdf} \tab  -\cr
 \tab \code{fromChipType} \tab  -\cr
 \tab \code{fromDataSet} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getData} \tab  -\cr
 \tab \code{getFields} \tab  -\cr
 \tab \code{getFragmentLengths} \tab  -\cr
 \tab \code{getFragmentStarts} \tab  -\cr
 \tab \code{getFragmentStops} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{readData} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
