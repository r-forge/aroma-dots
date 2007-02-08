%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QualityAssessmentSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{QualityAssessmentSet}
\docType{class}
\alias{QualityAssessmentSet}


\title{The QualityAssessmentSet class}

\description{
 Package:  aroma.affymetrix \cr
\bold{Class QualityAssessmentSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixCelSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{QualityAssessmentSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{QualityAssessmentSet}\cr
extends \link[aroma.affymetrix]{AffymetrixCelSet}\cr



}

\usage{QualityAssessmentSet(files=NULL, tags="*", ...)}

\arguments{
  \item{files}{A \code{\link[base]{list}} of \code{\link{QualityAssessmentFile}}:s.}
  \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
  \item{...}{Not used.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from AffymetrixCelSet}:\cr
[, [[, append, as, as.AffymetrixCelSet, averageQuantile, cna, extractAffyBatch, extractSnpFeatureSet, findByName, fromFiles, fromName, gcrmaSummary, getAverage, getAverageAsinh, getAverageFile, getAverageLog, getCdf, getData, getFullName, getIntensities, getTimestamps, getUnitIntensities, isDuplicated, nbrOfArrays, normalizeQuantile, plotDensity, range, readUnits, rmaSummary, setCdf

 \bold{Methods inherited from AffymetrixFileSet}:\cr
append, appendFiles, as, as.AffymetrixFileSet, as.list, clearCache, extract, fromFiles, getAlias, getFile, getFileSize, getFullName, getName, getNames, getPath, getPathnames, getTags, hasTag, hasTags, lapply, nbrOfFiles, sapply, seq, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}

\author{
  Ken Simpson (ksimpson[at]wehi.edu.au).
}
\keyword{classes}
