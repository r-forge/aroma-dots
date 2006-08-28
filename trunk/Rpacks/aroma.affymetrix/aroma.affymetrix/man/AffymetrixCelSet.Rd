%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixCelSet}
\docType{class}
\alias{AffymetrixCelSet}


 \title{The AffymetrixCelSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixCelSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFileSet}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{AffymetrixCelSet}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{ChipEffectSet}\cr

public static class \bold{AffymetrixCelSet}\cr
extends \link[aroma.affymetrix]{AffymetrixFileSet}\cr



  An AffymetrixCelSet object represents a data set of Affymetrix CEL files
  with \emph{identical} chip types.
 }

 \usage{AffymetrixCelSet(files=NULL, ...)}

 \arguments{
   \item{files}{A \code{\link[base]{list}} of \code{\link{AffymetrixCelFile}}:s.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{[} \tab  -\cr
 \tab \code{[[} \tab  -\cr
 \tab \code{applyToUnitIntensities} \tab  -\cr
 \tab \code{as} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:as.AffymetrixCelSet.AffymetrixCelSet]{as.AffymetrixCelSet}} \tab Coerce an object to an AffymetrixCelSet object.\cr
 \tab \code{\link[aroma.affymetrix:as.character.AffymetrixCelSet]{as.character}} \tab Returns a short string describing the Affymetrix data set.\cr
 \tab \code{clone} \tab  -\cr
 \tab \code{fromFiles} \tab  -\cr
 \tab \code{getAverageFile} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getCdf.AffymetrixCelSet]{getCdf}} \tab Gets the CDF structure for this CEL set.\cr
 \tab \code{\link[aroma.affymetrix:getIntensities.AffymetrixCelSet]{getIntensities}} \tab Gets cell intensities from a set of cells and a set of arrays.\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getUnitIntensities.AffymetrixCelSet]{getUnitIntensities}} \tab Gets cell signals for a subset of units and a subset of arrays.\cr
 \tab \code{\link[aroma.affymetrix:nbrOfArrays.AffymetrixCelSet]{nbrOfArrays}} \tab Gets the number of arrays in the file set.\cr
 \tab \code{readUnits} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:setCdf.AffymetrixCelSet]{setCdf}} \tab Sets the CDF structure for this CEL set.\cr
}


\bold{Methods inherited from AffymetrixFileSet}:\cr
as, as.AffymetrixFileSet, as.character, as.list, clearCache, clone, extract, fromFiles, getFile, getFileSize, getName, getNames, getPath, lapply, length, nbrOfFiles, seq

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
ds <- AffymetrixDataSet$fromFiles(path=path)

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
   \code{\link{AffymetrixCelFile}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
