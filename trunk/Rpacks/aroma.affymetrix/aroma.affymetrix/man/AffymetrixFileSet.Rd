%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixFileSet}
\docType{class}
\alias{AffymetrixFileSet}


 \title{The AffymetrixFileSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixFileSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{GenericDataFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaMicroarrayDataSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{AffymetrixFileSet}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffymetrixCelSet}, \link[aroma.affymetrix]{ChipEffectSet}, \link[aroma.affymetrix]{CnChipEffectSet}, \link[aroma.affymetrix]{DChipDcpSet}, \link[aroma.affymetrix]{ExonChipEffectSet}, \link[aroma.affymetrix]{FirmaSet}, \link[aroma.affymetrix]{ParameterCelSet}, \link[aroma.affymetrix]{QualityAssessmentSet}, \link[aroma.affymetrix]{ResidualSet}, \link[aroma.affymetrix]{SnpChipEffectSet}, \link[aroma.affymetrix]{WeightsSet}\cr

public static class \bold{AffymetrixFileSet}\cr
extends \link[aroma.core]{AromaMicroarrayDataSet}\cr



  An AffymetrixFileSet object represents a set of \code{\link{AffymetrixFile}}s
  with \emph{identical} chip types.
 }

 \usage{AffymetrixFileSet(files=NULL, ...)}

 \arguments{
   \item{files}{A \code{\link[base]{list}} of \code{\link{AffymetrixFile}}:s.}
   \item{...}{Arguments passed to \code{\link{GenericDataFileSet}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as} \tab  -\cr
 \tab \code{as.AffymetrixFileSet} \tab  -\cr
 \tab \code{fromFiles} \tab  -\cr
}


 \bold{Methods inherited from AromaMicroarrayDataSet}:\cr
getChipType, getPlatform, setAttributesBy, setAttributesBySampleAnnotationFile, setAttributesBySampleAnnotationSet, validate

 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, as.list, byName, clearCache, copyTo, equals, extract, findByName, fromFiles, getAlias, getFile, getFileClass, getFileSize, getFullName, getFullNames, getName, getNames, getPath, getPathnames, getTags, hasFile, hasTag, hasTags, indexOf, lapply, nbrOfFiles, sapply, seq, setAlias, setTags, validate

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
