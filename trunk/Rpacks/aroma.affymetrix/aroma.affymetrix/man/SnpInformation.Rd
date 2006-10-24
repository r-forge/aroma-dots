%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpInformation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{SnpInformation}
\docType{class}
\alias{SnpInformation}


 \title{The SnpInformation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class SnpInformation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{SnpInformation}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{DChipSnpInformation}\cr

public abstract static class \bold{SnpInformation}\cr
extends \link[aroma.affymetrix]{AffymetrixFile}\cr


 }

 \usage{SnpInformation(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:fromCdf.SnpInformation]{fromCdf}} \tab Static method to define a genome information set from a CDF.\cr
 \tab \code{\link[aroma.affymetrix:fromChipType.SnpInformation]{fromChipType}} \tab Static method to define a genome information set by chip type.\cr
 \tab \code{fromDataSet} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChipType.SnpInformation]{getChipType}} \tab Gets the chip type of this genome information set.\cr
 \tab \code{\link[aroma.affymetrix:getData.SnpInformation]{getData}} \tab Gets all or a subset of the genome information data.\cr
 \tab \code{getFields} \tab  -\cr
 \tab \code{getFragmentLengths} \tab  -\cr
 \tab \code{getFragmentStarts} \tab  -\cr
 \tab \code{getFragmentStops} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{readData} \tab  -\cr
}


\bold{Methods inherited from AffymetrixFile}:\cr
as.character, clone, copyFile, fromFile, getFilename, getFileSize, getFileType, getLabel, getName, getPath, getPathname, setLabel, stextLabel, stextLabels, stextSize

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson \url{http://www.braju.com/R/}}
\keyword{classes}
