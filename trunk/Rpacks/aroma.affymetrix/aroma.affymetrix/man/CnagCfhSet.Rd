%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnagCfhSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CnagCfhSet}
\docType{class}
\alias{CnagCfhSet}


 \title{The CnagCfhSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CnagCfhSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{GenericDataFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{CnagCfhSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CnagCfhSet}\cr
extends \link[aroma.affymetrix]{GenericDataFileSet}\cr



  An CnagCfhSet object represents a set of CNAG CFH files
  with \emph{identical} chip types.
 }

 \usage{CnagCfhSet(files=NULL, ...)}

 \arguments{
   \item{files}{A \code{\link[base]{list}} of \code{\link{CnagCfhFile}}:s.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{[} \tab  -\cr
 \tab \code{[[} \tab  -\cr
 \tab \code{append} \tab  -\cr
 \tab \code{as} \tab  -\cr
 \tab \code{as.CnagCfhSet} \tab  -\cr
 \tab \code{findByName} \tab  -\cr
 \tab \code{fromFiles} \tab  -\cr
 \tab \code{fromName} \tab  -\cr
 \tab \code{getAverage} \tab  -\cr
 \tab \code{getAverageAsinh} \tab  -\cr
 \tab \code{getAverageFile} \tab  -\cr
 \tab \code{getAverageLog} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getData} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getTimestamps} \tab  -\cr
 \tab \code{isDuplicated} \tab  -\cr
 \tab \code{nbrOfArrays} \tab  -\cr
 \tab \code{range} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
 \tab \code{setCdf} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, as, as.GenericDataFileSet, as.list, clearCache, copyTo, extract, fromFiles, getAlias, getFile, getFileSize, getFullName, getFullNames, getName, getNames, getPath, getPathnames, getTags, hasTag, hasTags, indexOf, lapply, nbrOfFiles, sapply, seq, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \seealso{
   \code{\link{CnagCfhFile}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
