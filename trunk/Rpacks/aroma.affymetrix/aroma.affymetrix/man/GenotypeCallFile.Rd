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
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{GenotypeCallFile}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{GenotypeCallXdrFile}\cr

public static class \bold{GenotypeCallFile}\cr
extends \link[aroma.affymetrix]{AffymetrixFile}\cr



  The abstract GenotypeCallFile class represents a file containing genotype
  calls for a given genotype parameter estimates.
 }

 \usage{GenotypeCallFile(..., cdf=NULL)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:000-LB-.GenotypeCallFile]{[}} \tab Gets the genotype calls.\cr
 \tab \code{\link[aroma.affymetrix:000-LB-LT_-.GenotypeCallFile]{[<-}} \tab Sets the genotype calls.\cr
 \tab \code{\link[aroma.affymetrix:fromFile.GenotypeCallFile]{fromFile}} \tab Static method defining a genotype call file object.\cr
 \tab \code{\link[aroma.affymetrix:getCdf.GenotypeCallFile]{getCdf}} \tab Gets the CDF.\cr
 \tab \code{\link[aroma.affymetrix:readUnits.GenotypeCallFile]{readUnits}} \tab Gets the genotype calls.\cr
 \tab \code{setCdf} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:updateUnits.GenotypeCallFile]{updateUnits}} \tab Updates the genotype calls.\cr
}


 \bold{Methods inherited from AffymetrixFile}:\cr
compareChecksum, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, readChecksum, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, testAttributes, validateChecksum, writeChecksum

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
