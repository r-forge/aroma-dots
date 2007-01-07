%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallXdrFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GenotypeCallXdrFile}
\docType{class}
\alias{GenotypeCallXdrFile}


 \title{The GenotypeCallXdrFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GenotypeCallXdrFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{GenotypeCallFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{GenotypeCallXdrFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{GenotypeCallXdrFile}\cr
extends \link[aroma.affymetrix]{GenotypeCallFile}\cr



  The abstract GenotypeCallXdrFile class represents a file containing genotype
  calls for a given genotype parameter estimates.
 }

 \usage{GenotypeCallXdrFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fromFile} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
 \tab \code{updateUnits} \tab  -\cr
}


 \bold{Methods inherited from GenotypeCallFile}:\cr
[, [<-, getCdf, getChipType, readUnits, updateUnits

 \bold{Methods inherited from AffymetrixFile}:\cr
fromFile, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\keyword{classes}
\keyword{internal}
