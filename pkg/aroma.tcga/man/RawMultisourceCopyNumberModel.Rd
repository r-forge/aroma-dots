%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RawMultisourceCopyNumberModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{RawMultisourceCopyNumberModel}
\docType{class}
\alias{RawMultisourceCopyNumberModel}


 \title{The RawMultisourceCopyNumberModel class}

 \description{
  Package:  aroma.tcga \cr
\bold{Class RawMultisourceCopyNumberModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.tcga]{MultisourceChromosomalModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.tcga]{MultisourceCopyNumberChromosomalModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{RawMultisourceCopyNumberModel}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{RawMultisourceCopyNumberModel}\cr
extends \emph{\link[aroma.tcga]{MultisourceCopyNumberChromosomalModel}}\cr



  This class represents an idenity multi-source copy-number model which
  returns the input as is.
 }

 \usage{RawMultisourceCopyNumberModel(...)}

 \arguments{
   \item{...}{Passed to the constructor of the superclass.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getSetTag} \tab  -\cr
}


 \bold{Methods inherited from MultisourceCopyNumberChromosomalModel}:\cr
extractRawCopyNumbers, fit

 \bold{Methods inherited from MultisourceChromosomalModel}:\cr
as.character, clearCache, extractRawGenomicSignals, fit, getAsteriskTags, getChipType, getChipTypes, getChromosomes, getDataSetList, getFullName, getFullNames, getGenome, getGenomeData, getGenomeFile, getName, getNames, getParentPath, getPath, getReportPath, getRootPath, getTags, indexOf, nbrOfArrays, nbrOfSources, setGenome, setName

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
