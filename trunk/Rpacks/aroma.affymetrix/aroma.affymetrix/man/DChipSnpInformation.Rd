%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipSnpInformation.R
% 
% on Sat Jan 06 01:51:18 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{DChipSnpInformation}
\docType{class}
\alias{DChipSnpInformation}


 \title{The DChipSnpInformation class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class DChipSnpInformation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{SnpInformation}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{DChipSnpInformation}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{DChipSnpInformation}\cr
extends \emph{\link[aroma.affymetrix]{SnpInformation}}\cr


 }

 \usage{DChipSnpInformation(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{SnpInformation}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fromChipType} \tab  -\cr
 \tab \code{readData} \tab  -\cr
}


\bold{Methods inherited from SnpInformation}:\cr
fromCdf, fromChipType, fromDataSet, getChipType, getData, getFields, getFragmentLengths, getFragmentStarts, getFragmentStops, nbrOfUnits, readData, readTableInternal, verify

\bold{Methods inherited from AffymetrixFile}:\cr
as.character, clone, copyFile, fromFile, getFilename, getFileSize, getFileType, getFullName, getLabel, getName, getPath, getPathname, getTags, setLabel, stextLabel, stextLabels, stextSize

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getInternalAddress, getStaticInstance, hasField, hashCode, isReferable, ll, load, names, newInstance, novirtual, objectSize, print, save, staticCode


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
