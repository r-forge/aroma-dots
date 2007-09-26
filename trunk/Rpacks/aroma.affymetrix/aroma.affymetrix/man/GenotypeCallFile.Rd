%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GenotypeCallFile}
\docType{class}
\alias{GenotypeCallFile}


 \title{The GenotypeCallFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GenotypeCallFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{GenotypeCallFile}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{GenotypeCallXdrFile}\cr

public static class \bold{GenotypeCallFile}\cr
extends \link[aroma.affymetrix]{GenericDataFile}\cr



  The abstract GenotypeCallFile class represents a file containing genotype
  calls for a given genotype parameter estimates.
 }

 \usage{GenotypeCallFile(..., cdf=NULL)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{GenericDataFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{[} \tab  -\cr
 \tab \code{[<-} \tab  -\cr
 \tab \code{fromFile} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
 \tab \code{setCdf} \tab  -\cr
 \tab \code{updateUnits} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{File format}{
  A GenotypeCallFile object is stored to file as a
  \code{\link[R.huge]{FileByteVector}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{classes}
\keyword{internal}
