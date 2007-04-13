%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CnagCfhFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CnagCfhFile}
\docType{class}
\alias{CnagCfhFile}


 \title{The CnagCfhFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class CnagCfhFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{CnagCfhFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CnagCfhFile}\cr
extends \link[aroma.affymetrix]{AffymetrixFile}\cr



  A CnagCfhFile object represents a single CNAG CFH file.
 }

 \usage{CnagCfhFile(..., cdf=NULL)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixFile}}.}
   \item{cdf}{An optional \code{\link{AffymetrixCdfFile}}}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{[} \tab  -\cr
 \tab \code{[[} \tab  -\cr
 \tab \code{clone} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fromFile.CnagCfhFile]{fromFile}} \tab Defines an CnagCfhFile object from a CFH file.\cr
 \tab \code{\link[aroma.affymetrix:getCdf.CnagCfhFile]{getCdf}} \tab Gets the CDF structure for this CEL file.\cr
 \tab \code{nbrOfCells} \tab  -\cr
 \tab \code{nbrOfSnps} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:setCdf.CnagCfhFile]{setCdf}} \tab Sets the CDF structure for this CEL file.\cr
}


 \bold{Methods inherited from AffymetrixFile}:\cr
compareChecksum, copyTo, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, readChecksum, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically part of an \code{\link{CnagCfhSet}}.
 }
\keyword{classes}
