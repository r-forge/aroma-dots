%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RmaProbeAffinityFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{RmaProbeAffinityFile}
\docType{class}
\alias{RmaProbeAffinityFile}


 \title{The RmaProbeAffinityFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class RmaProbeAffinityFile}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{RmaProbeAffinityFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{RmaProbeAffinityFile}\cr
extends \emph{\link[aroma.affymetrix]{ProbeAffinityFile}}\cr



  This class represents estimates of probe affinities in the
  RMA model, see \code{\link{AffymetrixRmaModel}}.
 }

 \usage{RmaProbeAffinityFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{ParameterCelFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


\bold{Methods inherited from ProbeAffinityFile}:\cr
createFrom, findUnitsTodo, getFirstCellIndices, readUnits, updateUnits, writeSpatial

\bold{Methods inherited from ParameterCelFile}:\cr
createFrom, decode, decodeUnit, decodeUnitGroup, encode, encodeUnit, encodeUnitGroup, readUnits, updateUnits

\bold{Methods inherited from AffymetrixCelFile}:\cr
[, [[, fromFile, getCdf, getChipType, getFields, getHeader, getRectangle, image270, normalizeQuantiles, readUnits, setCdf, transformAffine, updateUnits, writeSpatial

\bold{Methods inherited from AffymetrixFile}:\cr
as.character, fromFile, getFilename, getFileSize, getFileType, getName, getPath, getPathname

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
   This class subclasses the \code{\link{ProbeAffinityFile}} class so that the
   above estimates can be stored in a CEL file.  Please see that class for
   more details on the general idea behind encoding/decoding data in
   CEL files.

   The mapping/encoding is as follows:
   \itemize{
     \item{intensities}{Hold \code{phi} (probe affinities).}
     \item{stdvs}{Hold \code{stdvs} (SD of the probe affinities).}
     \item{pixels}{Currently just zeros.}
   }
   For more details about the file format, see also the source code of
   this class.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
