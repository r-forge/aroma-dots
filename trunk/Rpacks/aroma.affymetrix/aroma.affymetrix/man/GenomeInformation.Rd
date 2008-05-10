%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenomeInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GenomeInformation}
\docType{class}
\alias{GenomeInformation}


 \title{The GenomeInformation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GenomeInformation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{GenomeInformation}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffymetrixCsvGenomeInformation}, \link[aroma.affymetrix]{DChipGenomeInformation}, \link[aroma.affymetrix]{UgpGenomeInformation}\cr

public abstract static class \bold{GenomeInformation}\cr
extends \link[aroma.core]{GenericDataFile}\cr


 }

 \usage{GenomeInformation(..., .verify=TRUE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{GenericDataFile}}.}
   \item{.verify}{For internal use only.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byChipType} \tab  -\cr
 \tab \code{fromCdf} \tab  -\cr
 \tab \code{fromChipType} \tab  -\cr
 \tab \code{fromDataSet} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getChromosomes} \tab  -\cr
 \tab \code{getChromosomeStats} \tab  -\cr
 \tab \code{getData} \tab  -\cr
 \tab \code{getPositions} \tab  -\cr
 \tab \code{getUnitIndices} \tab  -\cr
 \tab \code{getUnitsOnChromosome} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{plotDensity} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
