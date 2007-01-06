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
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{GenotypeCallSet}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{GenotypeCallXdrSet}\cr

public static class \bold{GenotypeCallSet}\cr
extends \link[aroma.affymetrix]{AffymetrixFileSet}\cr



  The GenotypeCallSet class represents a set of genotype-call files.
 }

 \usage{GenotypeCallSet(files=NULL, ...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixFileSet}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{[} \tab  -\cr
 \tab \code{[[} \tab  -\cr
 \tab \code{as} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:as.GenotypeCallSet.GenotypeCallSet]{as.GenotypeCallSet}} \tab Coerce an object to an GenotypeCallSet object.\cr
 \tab \code{fromFiles} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getCdf.GenotypeCallSet]{getCdf}} \tab Gets the CDF structure for this CEL set.\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:nbrOfArrays.GenotypeCallSet]{nbrOfArrays}} \tab Gets the number of arrays in the file set.\cr
 \tab \code{readUnits} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:setCdf.GenotypeCallSet]{setCdf}} \tab Sets the CDF structure for this genotype call set.\cr
}


 \bold{Methods inherited from AffymetrixFileSet}:\cr
append, appendFiles, as, as.AffymetrixFileSet, as.list, clearCache, extract, fromFiles, getFile, getFileSize, getFullName, getName, getNames, getPath, getPathnames, getTags, lapply, nbrOfFiles, seq, setName, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\keyword{classes}
\keyword{internal}
