%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenericDataFileSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GenericDataFileSet}
\docType{class}
\alias{GenericDataFileSet}


 \title{The GenericDataFileSet class}

 \description{
  Package:  aroma.core \cr
\bold{Class GenericDataFileSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{GenericDataFileSet}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.core]{AromaMicroarrayDataSet}, \link[aroma.core]{AromaTabularBinarySet}, \link[aroma.core]{GenericTabularFileSet}, \link[aroma.core]{SampleAnnotationSet}, \link[aroma.core]{TabularTextFileSet}\cr

public static class \bold{GenericDataFileSet}\cr
extends \link[R.oo]{Object}\cr



  An GenericDataFileSet object represents a set of \code{\link{GenericDataFile}}s
  with \emph{identical} chip types.
 }

 \usage{GenericDataFileSet(files=NULL, tags="*", alias=NULL, ...)}

 \arguments{
   \item{files}{A \code{\link[base]{list}} of \code{\link{GenericDataFile}}:s.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be used for this file set.
      The string \code{"*"} indicates that it should be replaced by the
      tags part of the file set pathname.}
   \item{alias}{A \code{\link[base]{character}} string specifying a name alias overriding the
      name inferred from the pathname.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{append} \tab  -\cr
 \tab \code{appendFiles} \tab  -\cr
 \tab \code{as.list} \tab  -\cr
 \tab \code{byName} \tab  -\cr
 \tab \code{clearCache} \tab  -\cr
 \tab \code{equals} \tab  -\cr
 \tab \code{extract} \tab  -\cr
 \tab \code{findByName} \tab  -\cr
 \tab \code{fromFiles} \tab  -\cr
 \tab \code{getAlias} \tab  -\cr
 \tab \code{getFile} \tab  -\cr
 \tab \code{getFileClass} \tab  -\cr
 \tab \code{getFileSize} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getFullNames} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getNames} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getPathnames} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{hasFile} \tab  -\cr
 \tab \code{hasTag} \tab  -\cr
 \tab \code{hasTags} \tab  -\cr
 \tab \code{indexOf} \tab  -\cr
 \tab \code{lapply} \tab  -\cr
 \tab \code{nbrOfFiles} \tab  -\cr
 \tab \code{sapply} \tab  -\cr
 \tab \code{seq} \tab  -\cr
 \tab \code{setAlias} \tab  -\cr
 \tab \code{setTags} \tab  -\cr
 \tab \code{validate} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
