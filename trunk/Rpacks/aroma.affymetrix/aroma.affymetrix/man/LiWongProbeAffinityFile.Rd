%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  LiWongProbeAffinityFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{LiWongProbeAffinityFile}
\docType{class}
\alias{LiWongProbeAffinityFile}


 \title{The LiWongProbeAffinityFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class LiWongProbeAffinityFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ParameterCelFile}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeAffinityFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{LiWongProbeAffinityFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{LiWongProbeAffinityFile}\cr
extends \emph{\link[aroma.affymetrix]{ProbeAffinityFile}}\cr



  This class represents estimates of probe affinities in the
  Li \& Wong (2001) model.
 }

 \usage{LiWongProbeAffinityFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ParameterCelFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


\bold{Methods inherited from ProbeAffinityFile}:\cr
findUnitsTodo, getFirstCellIndices, readUnits, updateUnits, writeSpatial

\bold{Methods inherited from ParameterCelFile}:\cr
decode, decodeUnit, decodeUnitGroup, encode, encodeUnit, encodeUnitGroup, readUnits, updateUnits

\bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, annotateMvsA, calcMvsA, clearData, clone, createFrom, fromFile, getCdf, getChipType, getData, getHeader, getRectangle, image270, nbrOfCells, normalizeQuantiles, plotMvsA, plotMvsPosition, plotMvsX, readUnits, setCdf, smoothScatterMvsA, stextChipType, transformAffine, updateUnits, writeSpatial

\bold{Methods inherited from AffymetrixFile}:\cr
as.character, clone, copyFile, fromFile, getFilename, getFileSize, getFileType, getLabel, getName, getPath, getPathname, setLabel, stextLabel, stextLabels

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Stored parameters}{
   For any given unit group consisting of L probes (L >= 2), the
   following information is stored:
   \itemize{
     \item{phi}{(L \code{\link[base]{double}}s [floats] in (0,Inf)) The probe affinities,
       one for each probe.}
     \item{stdvs}{(L \code{\link[base]{double}}s [floats] in (0,Inf)) The standard deviation
       of the probe-affinity estimates, one for each probe.}
     \item{outliers}{(L logicals) Specifies if a probe was considered to
       be an outlier or not, one for each probe.}
     \item{nbrOfIterations}{(an \code{\link[base]{integer}} in [0,Inf]) Number of iterations
       before the algorithm converged.
       A \code{0} indicates that this unit group has not been fitted.}
     \item{converged}{(a \code{\link[base]{logical}}) If \code{\link[base:logical]{FALSE}}, the algorithm did not
       converge when fitting the phis (and thetas).}
     \item{convergedOutliers}{(a \code{\link[base]{logical}}) If \code{\link[base:logical]{FALSE}}, the algorithm
       did not converge in deciding what are outliers.}
   }

   For more details about these parameters, see \code{\link[affy]{fit.li.wong}}.
 }

 \section{File format}{
   This class subclasses \code{\link{ParameterCelFile}} so that the above
   estimates can be stored in a CEL file.  Please see that class for
   more details on the general idea behind encoding/decoding data in
   CEL files.

   The mapping/encoding is as follows:
   \itemize{
     \item{intensities}{Hold \code{phi} (probe affinities).}
     \item{stdvs}{Hold \code{stdvs} (SD of the probe affinities).}
     \item{pixels}{Hold \code{outliers}, \code{nbrOfIterations},
       \code{converged}, and \code{convergedOutliers}.
       This is done by letting the sign code for \code{outliers} such that
       a negative sign corresponds to \code{\link[base:logical]{TRUE}} and a positive to \code{\link[base:logical]{FALSE}}.
       Then the absolute value of \code{pixels[1]} equals
       \code{nbrOfIterations}.
       The \code{converged} is encoded as the second bit (module 2) of
       \code{pixel[2]}, and \code{convergedOutliers} as the third
       bit (modulo 4).  Note that we know that there are at least two
       probes, because otherwise we cannot fit the model.}
   }
   For more details about the file format, see also the source code of
   this class.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically obtain through the
   \code{getProbeAffinities()} method for the
   \code{\link{AffymetrixLiWongModel}} class.
 }
\keyword{classes}
