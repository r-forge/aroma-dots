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
\code{~~+--}\code{AffymetrixFileSet}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffymetrixCelSet}, \link[aroma.affymetrix]{ChipEffectSet}, \link[aroma.affymetrix]{CnChipEffectSet}, \link[aroma.affymetrix]{GdasAnnotationSet}, \link[aroma.affymetrix]{GenotypeCallSet}, \link[aroma.affymetrix]{GenotypeCallXdrSet}, \link[aroma.affymetrix]{SnpChipEffectSet}\cr

public static class \bold{AffymetrixFileSet}\cr
extends \link[R.oo]{Object}\cr



  An AffymetrixFileSet object represents a set of \code{\link{AffymetrixFile}}s
  with \emph{identical} chip types.
 }

 \usage{AffymetrixFileSet(files=NULL, tags="*", ...)}

 \arguments{
   \item{files}{A \code{\link[base]{list}} of \code{\link{AffymetrixFile}}:s.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be used for this data set.
      The string \code{"*"} indicates that it should be replaced by the
      tags part of the data set pathname.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{append} \tab  -\cr
 \tab \code{appendFiles} \tab  -\cr
 \tab \code{as} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:as.AffymetrixFileSet.AffymetrixFileSet]{as.AffymetrixFileSet}} \tab Coerce an object to an AffymetrixFileSet object.\cr
 \tab \code{\link[aroma.affymetrix:as.character.AffymetrixFileSet]{as.character}} \tab Returns a short string describing the Affymetrix file set.\cr
 \tab \code{\link[aroma.affymetrix:as.list.AffymetrixFileSet]{as.list}} \tab Returns the files of the Affymetrix file set.\cr
 \tab \code{clearCache} \tab  -\cr
 \tab \code{clone} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:extract.AffymetrixFileSet]{extract}} \tab Extract a subset of the file set.\cr
 \tab \code{\link[aroma.affymetrix:fromFiles.AffymetrixFileSet]{fromFiles}} \tab Defines an AffymetrixFileSet object by searching for Affymetrix files.\cr
 \tab \code{getFile} \tab  -\cr
 \tab \code{getFiles} \tab  -\cr
 \tab \code{getFileSize} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getFullName.AffymetrixFileSet]{getFullName}} \tab Gets the full name of the file set.\cr
 \tab \code{\link[aroma.affymetrix:getName.AffymetrixFileSet]{getName}} \tab Gets the name of the file set.\cr
 \tab \code{\link[aroma.affymetrix:getNames.AffymetrixFileSet]{getNames}} \tab Gets the names of the files in the file set.\cr
 \tab \code{\link[aroma.affymetrix:getPath.AffymetrixFileSet]{getPath}} \tab Gets the path (directory) of the file set.\cr
 \tab \code{\link[aroma.affymetrix:getPathnames.AffymetrixFileSet]{getPathnames}} \tab Gets the pathnames of the files in the file set.\cr
 \tab \code{\link[aroma.affymetrix:getTags.AffymetrixFileSet]{getTags}} \tab Gets the tags of the file set.\cr
 \tab \code{\link[aroma.affymetrix:lapply.AffymetrixFileSet]{lapply}} \tab Applies a function to each of the Affymetrix file files.\cr
 \tab \code{\link[aroma.affymetrix:length.AffymetrixFileSet]{length}} \tab Gets the number of files in the file set.\cr
 \tab \code{nbrOfFiles} \tab  -\cr
 \tab \code{reorder} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:seq.AffymetrixFileSet]{seq}} \tab Gets an vector of file indices.\cr
 \tab \code{\link[aroma.affymetrix:setName.AffymetrixFileSet]{setName}} \tab Sets the name of the file set.\cr
 \tab \code{\link[aroma.affymetrix:setTags.AffymetrixFileSet]{setTags}} \tab Sets the tags of the file set.\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
