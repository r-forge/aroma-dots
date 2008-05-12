%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaFreqBCghBinarySet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AromaFreqBCghBinarySet}
\docType{class}
\alias{AromaFreqBCghBinarySet}


 \title{The AromaFreqBCghBinarySet class}

 \description{
  Package:  aroma.cn \cr
\bold{Class AromaFreqBCghBinarySet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{GenericDataFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTabularBinarySet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.cn]{AromaSignalBinarySet}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{AromaFreqBCghBinarySet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AromaFreqBCghBinarySet}\cr
extends \link[aroma.cn]{AromaSignalBinarySet}\cr



  An AromaFreqBCghBinarySet object represents a set of
  \code{\link{AromaFreqBCghBinaryFile}}s with \emph{identical} chip types.
 }

 \usage{AromaFreqBCghBinarySet(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AromaSignalBinarySet}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byName} \tab  -\cr
}


 \bold{Methods inherited from AromaSignalBinarySet}:\cr
byName, findByName, getChipType, getPlatform, validate

 \bold{Methods inherited from AromaTabularBinarySet}:\cr
extractMatrix

 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, as.list, byName, clearCache, copyTo, extract, findByName, fromFiles, getAlias, getFile, getFileClass, getFileSize, getFullName, getFullNames, getName, getNames, getPath, getPathnames, getTags, hasTag, hasTags, indexOf, lapply, nbrOfFiles, sapply, seq, setAlias, setTags, validate

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
