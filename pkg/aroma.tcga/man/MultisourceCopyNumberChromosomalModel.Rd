%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MultisourceCopyNumberChromosomalModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{MultisourceCopyNumberChromosomalModel}
\docType{class}
\alias{MultisourceCopyNumberChromosomalModel}


 \title{The MultisourceCopyNumberChromosomalModel class}

 \description{
  Package:  aroma.tcga \cr
\bold{Class MultisourceCopyNumberChromosomalModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.tcga]{MultisourceChromosomalModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{MultisourceCopyNumberChromosomalModel}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.tcga]{RawMultisourceCopyNumberModel}}\cr

public abstract static class \bold{MultisourceCopyNumberChromosomalModel}\cr
extends \emph{\link[aroma.tcga]{MultisourceChromosomalModel}}\cr



  This \emph{abstract} class represents a multi-source copy-number model.
 }

 \usage{MultisourceCopyNumberChromosomalModel(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{MultisourceChromosomalModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.tcga:fit.MultisourceCopyNumberChromosomalModel]{fit}} \tab Fits the model.\cr
}


 \bold{Methods inherited from MultisourceChromosomalModel}:\cr
as.character, clearCache, extractRawGenomicSignals, fit, getAsteriskTags, getChipType, getChipTypes, getChromosomes, getDataSetList, getFullName, getFullNames, getGenome, getGenomeData, getGenomeFile, getName, getNames, getParentPath, getPath, getReportPath, getRootPath, getTags, indexOf, nbrOfArrays, nbrOfSources, setGenome, setName

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
