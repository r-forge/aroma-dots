%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GenotypeCallSet}
\docType{class}
\alias{GenotypeCallSet}


 \title{The GenotypeCallSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GenotypeCallSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{GenericDataFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{GenotypeCallSet}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{GenotypeCallXdrSet}\cr

public static class \bold{GenotypeCallSet}\cr
extends \link[aroma.affymetrix]{GenericDataFileSet}\cr



  The GenotypeCallSet class represents a set of genotype-call files.
 }

 \usage{GenotypeCallSet(files=NULL, ...)}

 \arguments{
   \item{files}{A \code{\link[base]{list}} of \code{\link{GenotypeCallFile}}:s.}
   \item{...}{Arguments passed to \code{\link{GenericDataFileSet}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{[} \tab  -\cr
 \tab \code{[[} \tab  -\cr
 \tab \code{as} \tab  -\cr
 \tab \code{as.GenotypeCallSet} \tab  -\cr
 \tab \code{fromFiles} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{nbrOfArrays} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
 \tab \code{setCdf} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, as, as.GenericDataFileSet, as.list, clearCache, copyTo, extract, fromFiles, getAlias, getFile, getFileSize, getFullName, getFullNames, getName, getNames, getPath, getPathnames, getTags, hasTag, hasTags, indexOf, lapply, nbrOfFiles, sapply, seq, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\keyword{classes}
\keyword{internal}
