%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastaReferenceSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{FastaReferenceSet}
\docType{class}
\alias{FastaReferenceSet}


\title{The FastaReferenceSet class}

\description{
 Package:  aroma.seq \cr
\bold{Class FastaReferenceSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFileSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{FastaReferenceSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{FastaReferenceSet}\cr
extends \link[R.filesets]{GenericDataFileSet}\cr



 An FastaReferenceSet object represents a set of \code{\link{FastaReferenceFile}}:s.
}

\usage{
FastaReferenceSet(files=NULL, ...)
}

\arguments{
  \item{files}{A \code{\link[base]{list}} of \code{\link{FastaReferenceFile}}:s.}
  \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFileSet}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byPath} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFileSet}:\cr
append, appendFiles, appendFullNamesTranslator, appendFullNamesTranslatorByNULL, appendFullNamesTranslatorByTabularTextFile, appendFullNamesTranslatorByTabularTextFileSet, appendFullNamesTranslatorBydata.frame, appendFullNamesTranslatorByfunction, appendFullNamesTranslatorBylist, as.character, as.list, byName, byPath, clearCache, clearFullNamesTranslator, clone, copyTo, equals, extract, findByName, getChecksum, getChecksumObjects, getDefaultFullName, getFile, getFileClass, getFileSize, getFiles, getFullNames, getNames, getOneFile, getPath, getPathnames, getSubdirs, hasFile, indexOf, nbrOfFiles, resetFullNames, seq, setFullNamesTranslator, sortBy, update2, updateFullName, updateFullNames, validate, getFullNameTranslatorSet, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, findRspReportTemplate, asThis


}

\author{Taku Tokuyasu}
\keyword{classes}