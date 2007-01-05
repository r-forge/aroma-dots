%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GladModel.R
% 
% on Sat Jan 06 01:52:50 2007.
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



  This class represents the GLAD [1] model.
 }

 \usage{GladModel(ces=NULL, reference=NULL, tags="*", ...)}

 \arguments{
   \item{ces}{A \code{\link{CnChipEffectSet}}.}
   \item{reference}{A \code{\link{CnChipEffectFile}}.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fit} \tab  -\cr
 \tab \code{getChipEffects} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChromosomes.GladModel]{getChromosomes}} \tab Gets the chromosomes available.\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getGenomeInformation} \tab  -\cr
 \tab \code{getLog2Ratios} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getReference} \tab  -\cr
 \tab \code{getRegions} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{nbrOfArrays} \tab  -\cr
 \tab \code{plot} \tab  -\cr
 \tab \code{setReference} \tab  -\cr
 \tab \code{writeRegions} \tab  -\cr
 \tab \code{writeWig} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getInternalAddress, getStaticInstance, hasField, hashCode, isReferable, ll, load, names, newInstance, novirtual, objectSize, print, save, staticCode


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
  [1] Hupe P et al. \emph{Analysis of array CGH data: from signal ratio to
      gain and loss of DNA regions}. Bioinformatics, 2004, 20, 3413-3422.\cr
 }
\keyword{classes}
