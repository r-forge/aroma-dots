%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QuantileNormalization.R
% 
% on Sat Jan 06 01:56:11 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{QuantileNormalization}
\docType{class}
\alias{QuantileNormalization}


 \title{The QuantileNormalization class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class QuantileNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{Preprocessing}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbePreprocessing}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{QuantileNormalization}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{OligoQuantileNormalization}\cr

public static class \bold{QuantileNormalization}\cr
extends \emph{\link[aroma.affymetrix]{ProbePreprocessing}}\cr



  This class represents a normalization function that transforms the
  probe-level signals towards the same empirical distribution.
 }

 \usage{QuantileNormalization(..., subsetToUpdate=NULL, typesToUpdate=NULL, targetDistribution=NULL, subsetToAvg=subsetToUpdate, typesToAvg=typesToUpdate)}

 \arguments{
   \item{...}{Arguments passed to the constructor of \code{\link{ProbePreprocessing}}.}
   \item{subsetToUpdate}{The probes to be updated.
     If \code{\link[base]{NULL}}, all probes are updated.}
   \item{typesToUpdate}{Types of probes to be updated.}
   \item{targetDistribution}{A \code{\link[base]{numeric}} \code{\link[base]{vector}}.  The empirical
     distribution to which all arrays should be normalized to.}
   \item{subsetToAvg}{The probes to calculate average empirical
     distribution over.  If a single \code{\link[base]{numeric}} in (0,1), then this
     fraction of all probes will be used.
     If \code{\link[base]{NULL}}, all probes are considered.}
   \item{typesToAvg}{Types of probes to be used when calculating the
     average empirical distribution.
     If \code{"pm"} and \code{"mm"} only perfect-match and mismatch
     probes are used, respectively. If \code{"pmmm"} both types are used.
   }
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:process.QuantileNormalization]{process}} \tab Normalizes the data set.\cr
}


\bold{Methods inherited from ProbePreprocessing}:\cr
getRootPath

\bold{Methods inherited from Preprocessing}:\cr
as.character, getAsteriskTag, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputFiles, getParameters, getParametersAsString, getPath, getRootPath, getTags, isDone, process, setTags

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getInternalAddress, getStaticInstance, hasField, hashCode, isReferable, ll, load, names, newInstance, novirtual, objectSize, print, save, staticCode


 }

 \examples{\dontrun{
   for (zzz in 0) {

# Setup verbose output
verbose <- Arguments$getVerbose(-2)
timestampOn(verbose)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define an example dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Find any dataset
path <- findCelSet(minCount=5)
if (is.null(path))
  break

ds <- AffymetrixCelSet$fromFiles(path)
print(ds)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
norm <- QuantileNormalization(ds, subsetToAvg=1/3)
dsQN <- process(norm, verbose=verbose)
print(dsQN)

} # for (zzz in 0)
rm(zzz)

 }}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
