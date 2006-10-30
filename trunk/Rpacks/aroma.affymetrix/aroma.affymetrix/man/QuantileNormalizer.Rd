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

 \usage{QuantileNormalizer(dataSet=NULL, subsetToUpdate=NULL, typesToUpdate=NULL, targetDistribution=NULL, subsetToAvg=subsetToUpdate, typesToAvg=typesToUpdate, rootPath="normQuantile", ...)}

 \arguments{
   \item{dataSet}{A \code{\link{AffymetrixCelSet}}.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{getInputDataSet} \tab  -\cr
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{getOutputPath} \tab  -\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{getTargetDistribution} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:isDone.QuantileNormalizer]{isDone}} \tab Checks if the data set is normalized or not.\cr
 \tab \code{\link[aroma.affymetrix:process.QuantileNormalizer]{process}} \tab Normalizes the data set.\cr
 \tab \code{setRootPath} \tab  -\cr
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

 \author{Henrik Bengtsson \url{http://www.braju.com/R/}}
\keyword{classes}
