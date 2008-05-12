%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaMicroarrayDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AromaMicroarrayDataSet}
\docType{class}
\alias{AromaMicroarrayDataSet}


 \title{The AromaMicroarrayDataSet class}

 \description{
  Package:  aroma.core \cr
\bold{Class AromaMicroarrayDataSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{GenericDataFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{AromaMicroarrayDataSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AromaMicroarrayDataSet}\cr
extends \link[aroma.core]{GenericDataFileSet}\cr



  An AromaMicroarrayDataSet object represents a set of
  \code{\link{AromaMicroarrayDataFile}}s with \emph{identical} chip types.
 }

 \usage{AromaMicroarrayDataSet(files=NULL, ...)}

 \arguments{
   \item{files}{A \code{\link[base]{list}} of \code{\link{AromaMicroarrayDataFile}}:s.}
   \item{...}{Arguments passed to \code{\link{GenericDataFileSet}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getPlatform} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, as, as.GenericDataFileSet, as.list, byName, clearCache, copyTo, extract, findByName, fromFiles, getAlias, getFile, getFileClass, getFileSize, getFullName, getFullNames, getName, getNames, getPath, getPathnames, getTags, hasTag, hasTags, indexOf, lapply, nbrOfFiles, sapply, seq, setAlias, setTags, validate

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
