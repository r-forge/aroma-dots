%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QuantileNormalizer.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{QuantileNormalizer}
\docType{class}
\alias{QuantileNormalizer}


 \title{The QuantileNormalizer class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class QuantileNormalizer}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{QuantileNormalizer}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{QuantileNormalizer}\cr
extends \link[R.oo]{Object}\cr



  This class represents a normalization function that transforms the
  probe-level signals towards the same empirical distribution.
 }

 \usage{QuantileNormalizer(dataSet=NULL, subsetToUpdate=NULL, typesToUpdate=NULL, targetDistribution=NULL, subsetToAvg=subsetToUpdate, typesToAvg=typesToUpdate, tags="*", ...)}

 \arguments{
   \item{dataSet}{A \code{\link{AffymetrixCelSet}}.}
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
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be appended to the tags of
      the input data set.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getFullName.QuantileNormalizer]{getFullName}} \tab Gets the full name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getInputDataSet.QuantileNormalizer]{getInputDataSet}} \tab Gets the source data set.\cr
 \tab \code{\link[aroma.affymetrix:getName.QuantileNormalizer]{getName}} \tab Gets the name of the output data set.\cr
 \tab \code{\link[aroma.affymetrix:getOutputDataSet.QuantileNormalizer]{getOutputDataSet}} \tab Gets the normalized data set.\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getPath.QuantileNormalizer]{getPath}} \tab Gets the path of the output data set.\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getTags.QuantileNormalizer]{getTags}} \tab Gets the tags of the output data set.\cr
 \tab \code{getTargetDistribution} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:isDone.QuantileNormalizer]{isDone}} \tab Checks if the data set is normalized or not.\cr
 \tab \code{\link[aroma.affymetrix:process.QuantileNormalizer]{process}} \tab Normalizes the data set.\cr
 \tab \code{setTags} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


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
norm <- QuantileNormalizer(ds, subsetToAvg=1/3)
dsQN <- process(norm, verbose=verbose)
print(dsQN)

} # for (zzz in 0)
rm(zzz)

 }}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
