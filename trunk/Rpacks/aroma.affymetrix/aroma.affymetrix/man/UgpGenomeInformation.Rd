%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  UgpGenomeInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{UgpGenomeInformation}
\docType{class}
\alias{UgpGenomeInformation}


 \title{The UgpGenomeInformation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class UgpGenomeInformation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{GenomeInformation}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{UgpGenomeInformation}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{UgpGenomeInformation}\cr
extends \emph{\link[aroma.affymetrix]{GenomeInformation}}\cr



  This class represents Aroma UGP genome information files.
 }

 \usage{UgpGenomeInformation(..., .ugp=NULL, .verify=TRUE)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{GenomeInformation}}.}
   \item{.ugp}{For internal use only.}
   \item{.verify}{For internal use only.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{clearCache} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fromChipType.UgpGenomeInformation]{fromChipType}} \tab Defines a UgpGenomeInformation object by chip type.\cr
 \tab \code{getChromosomes} \tab  -\cr
 \tab \code{getData} \tab  -\cr
 \tab \code{getUnitsOnChromosome} \tab  -\cr
 \tab \code{readData} \tab  -\cr
}


 \bold{Methods inherited from GenomeInformation}:\cr
fromCdf, fromChipType, fromDataSet, getChipType, getChromosomes, getChromosomeStats, getData, getPositions, getUnitIndices, getUnitsOnChromosome, nbrOfUnits, plotDensity, readData

 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, isFile, readChecksum, renameTo, renameToUpperCaseExt, setAlias, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
