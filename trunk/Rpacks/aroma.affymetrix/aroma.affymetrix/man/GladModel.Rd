%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GladModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GladModel}
\docType{class}
\alias{GladModel}


 \title{The GladModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GladModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{GladModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{GladModel}\cr
extends \link[R.oo]{Object}\cr



  This class represents the Gain and Loss Analysis of DNA regions
  (GLAD) model [1].
  This class can model chip-effect estimates obtained from multiple
  chip types, and not all samples have to be available on all chip types.
 }

 \usage{GladModel(cesList=NULL, referenceList=NULL, tags="*", ...)}

 \arguments{
   \item{cesList}{A \code{\link[base]{list}} of \code{\link{ChipEffectSet}} objects.}
   \item{referenceList}{A \code{\link{ChipEffectFile}} objects.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fit} \tab  -\cr
 \tab \code{getArrays} \tab  -\cr
 \tab \code{getChipEffectFiles} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getChipTypes} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChromosomes.GladModel]{getChromosomes}} \tab Gets the chromosomes available.\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getListOfChipEffects} \tab  -\cr
 \tab \code{getListOfGenomeInformations} \tab  -\cr
 \tab \code{getListOfReferences} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getNames} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getReferenceFiles} \tab  -\cr
 \tab \code{getRegions} \tab  -\cr
 \tab \code{getTableOfArrays} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{nbrOfArrays} \tab  -\cr
 \tab \code{nbrOfChipTypes} \tab  -\cr
 \tab \code{plot} \tab  -\cr
 \tab \code{writeRegions} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \details{
   Data from multiple chip types are combined "as is".  This is based
   on the assumption that the relative chip effect estimates are
   non-biased (or at the equally biased across chip types).
   Note that in GLAD there is no way to down weight certain data points,
   which is why we can control for differences in variance across
   chip types.
 }

 \section{Benchmarking}{
   In high-density copy numbers analysis, the most time consuming step
   is fitting the GLAD model.  The complexity of the model grows
   more than linearly (squared? exponentially?) with the number of data
   points in the chromosome and sample being fitted.  This is why it
   take much more than twice the time to fit two chip types together
   than separately.
 }

 \section{Requirements}{
   This class requires genome information annotation files for
   every chip type.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
  \code{\link{GladModel}}.
 }

 \references{
  [1] Hupe P et al. \emph{Analysis of array CGH data: from signal ratio to
      gain and loss of DNA regions}. Bioinformatics, 2004, 20, 3413-3422.\cr
 }
\keyword{classes}
