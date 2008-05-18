%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenomeInformation.R
% 
% on Sat May 17 11:54:49 2008.
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
 \tab \code{\link[aroma.affymetrix:byChipType.GenomeInformation]{byChipType}} \tab Static method to define a genome information set by chip type.\cr
 \tab \code{\link[aroma.affymetrix:fromCdf.GenomeInformation]{fromCdf}} \tab Static method to define a genome information set from a CDF.\cr
 \tab \code{fromChipType} \tab  -\cr
 \tab \code{fromDataSet} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChipType.GenomeInformation]{getChipType}} \tab Gets the chip type of this genome information set.\cr
 \tab \code{getChromosomes} \tab  -\cr
 \tab \code{getChromosomeStats} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getData.GenomeInformation]{getData}} \tab Gets all or a subset of the genome information data.\cr
 \tab \code{\link[aroma.affymetrix:getPositions.GenomeInformation]{getPositions}} \tab Gets the physical positions for a set of units.\cr
 \tab \code{\link[aroma.affymetrix:getUnitIndices.GenomeInformation]{getUnitIndices}} \tab Gets unit indices ordered along the chromosome.\cr
 \tab \code{getUnitsOnChromosome} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:plotDensity.GenomeInformation]{plotDensity}} \tab Plots the density of SNPs for a given chromosome.\cr
}


 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getFullNameTranslator, getName, getPath, getPathname, getTags, gunzip, gzip, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setFullNameTranslator, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
