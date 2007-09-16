%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixFile}
\docType{class}
\alias{AffymetrixFile}


 \title{The abstract AffymetrixFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{GenericDataFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{AffymetrixFile}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffymetrixCdfFile}, \link[aroma.affymetrix]{AffymetrixCelFile}, \link[aroma.affymetrix]{AffymetrixProbeTabFile}, \link[aroma.affymetrix]{AffymetrixTsvFile}, \link[aroma.affymetrix]{ChipEffectFile}, \link[aroma.affymetrix]{CnChipEffectFile}, \link[aroma.affymetrix]{CnProbeAffinityFile}, \link[aroma.affymetrix]{ExonChipEffectFile}, \link[aroma.affymetrix]{ExonProbeAffinityFile}, \link[aroma.affymetrix]{FirmaFile}, \link[aroma.affymetrix]{GdasAnnotationFile}, \link[aroma.affymetrix]{ParameterCelFile}, \link[aroma.affymetrix]{ProbeAffinityFile}, \link[aroma.affymetrix]{QualityAssessmentFile}, \link[aroma.affymetrix]{ResidualFile}, \link[aroma.affymetrix]{SampleAnnotationFile}, \link[aroma.affymetrix]{SnpChipEffectFile}, \link[aroma.affymetrix]{SnpProbeAffinityFile}, \link[aroma.affymetrix]{WeightsFile}\cr

public static class \bold{AffymetrixFile}\cr
extends \link[aroma.affymetrix]{GenericDataFile}\cr



  An AffymetrixFile object represents a single Affymetrix file,
  e.g. an Affymetrix CEL file or an Affymetrix CDF file.
  Note that this class is abstract and can not be instanciated, but
  instead you have to use one of the subclasses or the generic
  \code{\link[aroma.affymetrix:fromFile.AffymetrixFile]{*fromFile}()} method.
 }

 \usage{AffymetrixFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{GenericDataFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from GenericDataFile}:\cr
compareChecksum, copyTo, equals, fromFile, getAlias, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, readChecksum, renameTo, renameToUpperCaseExt, setAlias, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically part of an \code{\link{AffymetrixFileSet}}.
 }
\keyword{classes}
