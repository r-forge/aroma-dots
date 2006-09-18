###########################################################################/**
# @RdocDocumentation "Non-documented objects"
#
# % The AffymetrixCelSet class
# @alias getSiblings
#
# % The AffymetrixDataFile class
# @alias addTransform
# @alias addTransform.AffymetrixDataFile
# @alias createOptimalReadMap
# @alias createOptimalReadMap.AffymetrixDataFile
# @alias fromFile
# @alias getApdMap
# @alias getApdMap.AffymetrixDataFile
# @alias getChipType
# @alias getChipType.AffymetrixApdFile
# @alias getProbes
# @alias getProbes.AffymetrixDataFile
# @alias hasApdMap
# @alias readUnits
#
# % The AffymetrixApdFile class
# @alias getFileType.AffymetrixApdFile
# @alias readUnits.AffymetrixApdFile
#
# % The AffymetrixCelFile class
# @alias getFileType.AffymetrixCelFile
#
# % The AffymetrixDataset class
# @alias as.AffymetrixDataset
# @alias compactApdMap
# @alias compactApdMap.AffymetrixDataset
# @alias createOptimalReadMap.AffymetrixDataset
# @alias getProbes.AffymetrixDataset
# @alias findAverageQuantile
# @alias findAverageQuantile.AffymetrixDataset
# @alias getApdMap.AffymetrixDataset
# @alias getFileType
# @alias getFileType.AffymetrixDataset
#
# @alias newInstance.Object
# @alias getSnpNames
#
# % The ApdMap class
# @alias as.character.ApdMap
# @alias equals.ApdMap
# @alias findApdMap.ApdMap
# @alias fromCdf
# @alias fromCdf.ApdMap
# @alias fromMapType
# @alias fromMapType.ApdMap
# @alias getChipType.ApdMap
# @alias getMapType.ApdMap
# @alias getReadMap
# @alias getReadMap.ApdMap
# @alias getWriteMap
# @alias getWriteMap.ApdMap
# @alias read
# @alias read.ApdMap
# @alias write
# @alias write.ApdMap
#
# % The AffymetrixModel class
# @alias as.character.AffymetrixModel
# @alias findUnitsDone
# @alias findUnitsDone.AffymetrixModel
# @alias fit
# @alias fit.AffymetrixModel
# @alias fromAffymetrixDataset
# @alias fromAffymetrixDataset.AffymetrixModel
# @alias getDataset
# @alias getDataset.AffymetrixModel
# @alias getFit
# @alias getFit.AffymetrixModel
# @alias getPath.AffymetrixModel
# @alias nbrOfArrays
# @alias nbrOfArrays.AffymetrixModel
# @alias nbrOfProbes.AffymetrixModel
# @alias nbrOfProbesets
# @alias nbrOfProbesets.AffymetrixModel
#
# @alias colMeans
# @alias colSums
# @alias doTransform
# @alias findCdf
# @alias findCdf.default
# @alias fitQuantileNormFcn
# @alias fromFiles
# @alias getFilename
# @alias getMapType
# @alias getMapType.AffymetrixApdFile
# @alias getMapType.AffymetrixCelFile
# @alias getNames
# @alias getProbeIndices
# @alias getProbeIntensities
# @alias getProbesetIntensities
# @alias getProbesetNames
# @alias getProbesets
# @alias getProbesets.AffymetrixDataFile
# @alias hasReadMap
# @alias hasTransforms
# @alias hasTransforms.AffymetrixDataFile
# @alias isTransforming
# @alias isTransforming.AffymetrixDataFile
# @alias lapply
# @alias lapply.default
# @alias nbrOfProbes
# @alias nbrOfProbes.AffymetrixApdFile
# @alias nbrOfProbes.AffymetrixCelFile
# @alias nbrOfSnps
# @alias readCdfUnitsWriteMap
# @alias readCdfUnitsWriteMap.default
# @alias readIntensities
# @alias readIntensities.AffymetrixApdFile
# @alias readIntensities.AffymetrixCelFile
# @alias setApdMap
# @alias setTransform
# @alias setTransform.AffymetrixDataFile
# @alias setTransformStatus
# @alias setTransformStatus.AffymetrixDataFile
# @alias stringTree
# @alias stringTree.character
# @alias transformOff
# @alias transformOff.AffymetrixDataFile
# @alias transformOn
# @alias transformOn.AffymetrixDataFile
# @alias transformProbeSignals
# @alias transformProbeSignals.AffymetrixDataFile
# @alias write.default
# @alias writeApd.AffymetrixDataFile
# % Calibration and normalization functions
#
# % Miscellaneous statistical functions
#
# % Other missing docs
# @eval "t <- readLines('999.missingdocs.txt'); t <- trim(unlist(strsplit(t, split=' '))); t <- t[nchar(t) > 0]; t2 <- gsub('\\[', '\\\\[', t); t <- unique(t); t <- sprintf('\\alias{%s}', t); paste(t, collapse='\n')"
#
# \description{
#   This page contains aliases for all "non-documented" objects that 
#   \code{R CMD check} detects in this package. 
#
#   Almost all of them are \emph{generic} functions that have specific 
#   document for the corresponding method coupled to a specific class. 
#   Other functions are re-defined by \code{setMethodS3()} to 
#   \emph{default} methods. Neither of these two classes are non-documented
#   in reality.
#   The rest are deprecated methods.
# }
#
# @author
#
# @keyword internal
#*/###########################################################################

############################################################################
# HISTORY:
# 2005-05-15
# o Created to please R CMD check.
############################################################################
