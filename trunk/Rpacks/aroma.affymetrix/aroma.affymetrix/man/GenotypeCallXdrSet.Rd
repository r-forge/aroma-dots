%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallXdrSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GenotypeCallXdrSet}
\docType{class}
\alias{GenotypeCallXdrSet}


 \title{The GenotypeCallXdrSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GenotypeCallXdrSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{GenotypeCallSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{GenotypeCallXdrSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{GenotypeCallXdrSet}\cr
extends \link[aroma.affymetrix]{GenotypeCallSet}\cr



  The GenotypeCallXdrSet class represents a GenotypeCallSet where the
  data files are stored in the binary XDR file format.

  \emph{You should consider this class temporary until another file format
  is defined for genotype calls. /HB 2006-12-14.}
 }

 \usage{GenotypeCallXdrSet(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{GenotypeCallSet}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fromFiles} \tab  -\cr
}


 \bold{Methods inherited from GenotypeCallSet}:\cr
[, [[, as, as.GenotypeCallSet, fromFiles, getCdf, getFullName, nbrOfArrays, readUnits, setCdf

 \bold{Methods inherited from AffymetrixFileSet}:\cr
append, appendFiles, as, as.AffymetrixFileSet, as.list, clearCache, extract, fromFiles, getAlias, getFile, getFileSize, getFullName, getName, getNames, getPath, getPathnames, getTags, hasTag, hasTags, lapply, nbrOfFiles, sapply, seq, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\keyword{classes}
\keyword{internal}
