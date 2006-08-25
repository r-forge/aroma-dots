%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixDataSet}
\docType{class}
\alias{AffymetrixDataSet}


 \title{The AffymetrixDataSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixDataSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{AffymetrixDataSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AffymetrixDataSet}\cr
extends \link[R.oo]{Object}\cr



  An AffymetrixDataSet object represents an Affymetrix data set consisting
  of one or several Affymetrix data files with \emph{identical} chip types.
  The purpose of this class is to keep memory usage down.
  For this reason, each time data is requested, internally it is read from
  file and optionally transformed via, say, normalization functions.
 }

 \usage{AffymetrixDataSet(dataFiles=NULL, ...)}

 \arguments{
   \item{dataFiles}{A \code{\link[base]{list}} of \code{\link{AffymetrixDataFile}}:s.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:000-LB-.AffymetrixDataSet]{[}} \tab Extract a subset of the data set.\cr
 \tab \code{\link[aroma.affymetrix:000-LB--LB-.AffymetrixDataSet]{[[}} \tab Gets a subset of the data files in the data set.\cr
 \tab \code{as} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:as.AffymetrixDataSet.AffymetrixDataSet]{as.AffymetrixDataSet}} \tab Coerce an object to an AffymetrixDataSet object.\cr
 \tab \code{\link[aroma.affymetrix:as.character.AffymetrixDataSet]{as.character}} \tab Returns a short string describing the Affymetrix data set.\cr
 \tab \code{\link[aroma.affymetrix:as.list.AffymetrixDataSet]{as.list}} \tab Returns the data files of the Affymetrix data set.\cr
 \tab \code{\link[aroma.affymetrix:averageQuantile.AffymetrixDataSet]{averageQuantile}} \tab Gets the average empirical distribution across all samples.\cr
 \tab \code{clearCache} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:extract.AffymetrixDataSet]{extract}} \tab Extract a subset of the data set.\cr
 \tab \code{\link[aroma.affymetrix:fromFiles.AffymetrixDataSet]{fromFiles}} \tab Defines an AffymetrixDataSet object by searching for Affymetrix files.\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getFileType} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getNames.AffymetrixDataSet]{getNames}} \tab Gets the names of the arrays in the data set.\cr
 \tab \code{\link[aroma.affymetrix:getPath.AffymetrixDataSet]{getPath}} \tab Gets the path (directory) of the data set.\cr
 \tab \code{\link[aroma.affymetrix:getProbeIntensities.AffymetrixDataSet]{getProbeIntensities}} \tab Gets probe intensities from a set of probes and a set of arrays.\cr
 \tab \code{\link[aroma.affymetrix:getProbesetIntensities.AffymetrixDataSet]{getProbesetIntensities}} \tab Gets probe signals for a subset of probesets and a subset of arrays.\cr
 \tab \code{\link[aroma.affymetrix:lapply.AffymetrixDataSet]{lapply}} \tab Applies a function to each of the Affymetrix data files.\cr
 \tab \code{\link[aroma.affymetrix:nbrOfArrays.AffymetrixDataSet]{nbrOfArrays}} \tab Gets the number of arrays in the data set.\cr
 \tab \code{\link[aroma.affymetrix:normalizeQuantile.AffymetrixDataSet]{normalizeQuantile}} \tab Normalizes samples to have the same empirical distribution.\cr
 \tab \code{\link[aroma.affymetrix:plotDensity.AffymetrixDataSet]{plotDensity}} \tab Plots the densities of all samples.\cr
 \tab \code{\link[aroma.affymetrix:seq.AffymetrixDataSet]{seq}} \tab Gets an vector of data file indices.\cr
 \tab \code{transformAffine} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:write.AffymetrixDataSet]{write}} \tab Writes the data files in this data set to file.\cr
 \tab \code{\link[aroma.affymetrix:writeApd.AffymetrixDataSet]{writeApd}} \tab Generates APD files for each of the data files in this data set.\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \examples{

for (zzz in 0) {

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Look for the 100K Hind data directory
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- filePath("cel/Xba/", expandLinks="any")
if (!isDirectory(path))
  break


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define a dataset object for all files therein
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds <- AffymetrixDataset$fromFiles(path=path)

# Keep at most sex arrays for this example
ds <- extract(ds, 1:min(6,nbrOfArrays(ds)))
print(ds)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsQN <- normalizeQuantile(ds, subsetToAvg=0.50,
                                  typesToUpdate="pm", verbose=TRUE)
print(dsQN)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot densities
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
layout(matrix(1:2, ncol=1))
xlim <- c(6,16)
ylim <- c(0,0.40)
plotDensity(ds, types="pm", xlim=xlim, ylim=ylim, main="Raw data")
plotDensity(dsQN, types="pm", xlim=xlim, ylim=ylim, main="Normalized data")

# Clean up
layout(1)
rm(path, xlim, ylim)
# rm(ds, dsQN)
} # for (zzz in 0)
}

 \seealso{
   \code{\link{AffymetrixApdFile}} and \code{\link{AffymetrixCelFile}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
