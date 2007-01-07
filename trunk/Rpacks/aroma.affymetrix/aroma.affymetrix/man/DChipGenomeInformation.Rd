%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipGenomeInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{DChipGenomeInformation}
\docType{class}
\alias{DChipGenomeInformation}


 \title{The DChipGenomeInformation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class DChipGenomeInformation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{GenomeInformation}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{DChipGenomeInformation}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{DChipGenomeInformation}\cr
extends \emph{\link[aroma.affymetrix]{GenomeInformation}}\cr


 }

 \usage{DChipGenomeInformation(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{GenomeInformation}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fromChipType} \tab  -\cr
 \tab \code{readData} \tab  -\cr
}


 \bold{Methods inherited from GenomeInformation}:\cr
fromCdf, fromChipType, fromDataSet, getChipType, getChromosomes, getChromosomeStats, getData, getPositions, getUnitIndices, getUnitsOnChromosome, nbrOfUnits, plotDensity, readData

 \bold{Methods inherited from AffymetrixFile}:\cr
fromFile, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
