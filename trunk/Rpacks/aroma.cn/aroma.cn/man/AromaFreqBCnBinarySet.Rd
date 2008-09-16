%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaFreqBCnBinarySet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AromaFreqBCnBinarySet}
\docType{class}
\alias{AromaFreqBCnBinarySet}


 \title{The AromaFreqBCnBinarySet class}

 \description{
  Package:  aroma.cn \cr
\bold{Class AromaFreqBCnBinarySet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{GenericDataFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{GenericTabularFileSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.core]{AromaTabularBinarySet}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.cn]{AromaSignalBinarySet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{AromaFreqBCnBinarySet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AromaFreqBCnBinarySet}\cr
extends \link[aroma.cn]{AromaSignalBinarySet}\cr



  An AromaFreqBCnBinarySet object represents a set of
  \code{\link{AromaFreqBCnBinaryFile}}s with \emph{identical} chip types.
 }

 \usage{AromaFreqBCnBinarySet(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AromaSignalBinarySet}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byName} \tab  -\cr
}


 \bold{Methods inherited from AromaSignalBinarySet}:\cr
byName, findByName, getAromaUgpFile, getChipType, getPlatform, validate

 \bold{Methods inherited from GenericTabularFileSet}:\cr
extractMatrix

 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, as.list, byName, clearCache, copyTo, equals, extract, findByName, fromFiles, getAlias, getFile, getFileClass, getFileSize, getFullName, getFullNames, getFullNameTranslator, getName, getNames, getPath, getPathnames, getTags, hasFile, hasTag, hasTags, indexOf, lapply, nbrOfFiles, sapply, seq, setAlias, setFullName, setFullNamesTranslator, setFullNameTranslator, setName, setTags, update2, validate

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
