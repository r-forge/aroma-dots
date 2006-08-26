%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeAffinityFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ProbeAffinityFile}
\docType{class}
\alias{ProbeAffinityFile}


 \title{The ProbeAffinityFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ProbeAffinityFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\emph{\code{ProbeAffinityFile}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{LiWongProbeAffinityFile}, \link[aroma.affymetrix]{RmaProbeAffinityFile}\cr

public abstract static class \bold{ProbeAffinityFile}\cr
extends \emph{\link[aroma.affymetrix]{ParameterCelFile}}\cr



  This class represents estimates of probe affinities in the
  RMA model.
 }

 \usage{ProbeAffinityFile(..., model=c("pm"))}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ParameterCelFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{createFrom} \tab  -\cr
 \tab \code{findUnitsTodo} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
}


\bold{Methods inherited from ParameterCelFile}:\cr
createFrom, decode, decodeUnit, decodeUnitGroup, encode, encodeUnit, encodeUnitGroup, readUnits, updateUnits

\bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, fromFile, getCdf, getChipType, getHeader, getIntensities, getRectangle, image270, readUnits, setCdf, updateUnits

\bold{Methods inherited from AffymetrixFile}:\cr
as.character, fromFile, getFilename, getFileSize, getFileType, getName, getPath, getPathname

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically obtained through the
   \code{getProbeAffinities()} method for the \code{\link{ProbeLevelModel}} class.
 }
\keyword{classes}
