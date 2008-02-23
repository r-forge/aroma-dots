%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  DChipDcpSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{DChipDcpSet}
\docType{class}
\alias{DChipDcpSet}


 \title{The DChipDcpSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class DChipDcpSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{GenericDataFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{DChipDcpSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{DChipDcpSet}\cr
extends \link[aroma.affymetrix]{AffymetrixFileSet}\cr



  A DChipDcpSet object represents a set of DChip DCP files
  for \emph{identical} chip types.
 }

 \usage{DChipDcpSet(files=NULL, ...)}

 \arguments{
   \item{files}{A \code{\link[base]{list}} of \code{\link{DChipDcpFile}}:s.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:as.DChipDcpSet.DChipDcpSet]{as.DChipDcpSet}} \tab Coerce an object to an DChipDcpSet object.\cr
 \tab \code{byName} \tab  -\cr
 \tab \code{findByName} \tab  -\cr
 \tab \code{fromFiles} \tab  -\cr
 \tab \code{fromName} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:nbrOfArrays.DChipDcpSet]{nbrOfArrays}} \tab Gets the number of arrays in the file set.\cr
}


 \bold{Methods inherited from AffymetrixFileSet}:\cr
as, as.AffymetrixFileSet, clearCache, fromFiles

 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, as, as.GenericDataFileSet, as.list, clearCache, copyTo, extract, fromFiles, getAlias, getFile, getFileSize, getFullName, getFullNames, getName, getNames, getPath, getPathnames, getTags, hasTag, hasTags, indexOf, lapply, nbrOfFiles, sapply, seq, setAlias, setTags

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \seealso{
   \code{\link{DChipDcpFile}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
