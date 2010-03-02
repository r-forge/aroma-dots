%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MultisourceChromosomalModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{MultisourceChromosomalModel}
\docType{class}
\alias{MultisourceChromosomalModel}


 \title{The MultisourceChromosomalModel class}

 \description{
  Package:  aroma.tcga \cr
\bold{Class MultisourceChromosomalModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{MultisourceChromosomalModel}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.tcga]{MultisourceCopyNumberChromosomalModel}}, \emph{\link[aroma.tcga]{RawMultisourceCopyNumberModel}}\cr

public abstract static class \bold{MultisourceChromosomalModel}\cr
extends \link[R.oo]{Object}\cr



  This \emph{abstract} class represents a chromosomal model based on
  data from multiple sources (platforms, labs, methods, ...).
 }

 \usage{MultisourceChromosomalModel(dsList=NULL, tags="*", genome="Human", ...)}

 \arguments{
   \item{dsList}{A \code{\link[base]{list}} of \code{\link{AromaUnitSignalBinarySet}}.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
   \item{genome}{A \code{\link[base]{character}} string specifying what genome is process.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{clearCache} \tab  -\cr
 \tab \code{fit} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getChipTypes} \tab  -\cr
 \tab \code{\link[aroma.tcga:getChromosomes.MultisourceChromosomalModel]{getChromosomes}} \tab Gets the chromosomes available.\cr
 \tab \code{getDataSetList} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getFullNames} \tab  -\cr
 \tab \code{getGenome} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{\link[aroma.tcga:getNames.MultisourceChromosomalModel]{getNames}} \tab Gets the names of the arrays.\cr
 \tab \code{getParentPath} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{indexOf} \tab  -\cr
 \tab \code{\link[aroma.tcga:nbrOfArrays.MultisourceChromosomalModel]{nbrOfArrays}} \tab Gets the number of arrays.\cr
 \tab \code{nbrOfSources} \tab  -\cr
 \tab \code{setGenome} \tab  -\cr
 \tab \code{setName} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \section{Requirements}{
   This class requires genome information annotation files for
   every chip type.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
