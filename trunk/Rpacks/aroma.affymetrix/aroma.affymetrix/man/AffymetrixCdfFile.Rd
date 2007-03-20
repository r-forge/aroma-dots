%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCdfFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixCdfFile}
\docType{class}
\alias{AffymetrixCdfFile}


 \title{The AffymetrixCdfFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixCdfFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{AffymetrixCdfFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AffymetrixCdfFile}\cr
extends \link[aroma.affymetrix]{AffymetrixFile}\cr



  An AffymetrixCdfFile object represents a generic Affymetrix CDF file.
 }

 \usage{AffymetrixCdfFile(...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:convert.AffymetrixCdfFile]{convert}} \tab Converts a CDF into the same CDF but with another format.\cr
 \tab \code{\link[aroma.affymetrix:fromChipType.AffymetrixCdfFile]{fromChipType}} \tab Defines an AffymetrixCdfFile object by chip type.\cr
 \tab \code{\link[aroma.affymetrix:fromFile.AffymetrixCdfFile]{fromFile}} \tab Defines an AffymetrixCdfFile object from a CDF file.\cr
 \tab \code{\link[aroma.affymetrix:getCellIndices.AffymetrixCdfFile]{getCellIndices}} \tab Gets the cell indices unit by unit.\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getDimension} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getGenomeInformation.AffymetrixCdfFile]{getGenomeInformation}} \tab Gets genome information for this chip type.\cr
 \tab \code{\link[aroma.affymetrix:getUnitNames.AffymetrixCdfFile]{getUnitNames}} \tab Gets the names of each unit.\cr
 \tab \code{\link[aroma.affymetrix:indexOf.AffymetrixCdfFile]{indexOf}} \tab Gets the indices of units by their names.\cr
 \tab \code{\link[aroma.affymetrix:isPm.AffymetrixCdfFile]{isPm}} \tab Checks which cells (probes) are PMs and not.\cr
 \tab \code{nbrOfCells} \tab  -\cr
 \tab \code{nbrOfColumns} \tab  -\cr
 \tab \code{nbrOfQcUnits} \tab  -\cr
 \tab \code{nbrOfRows} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:readUnits.AffymetrixCdfFile]{readUnits}} \tab Reads CDF data unit by unit.\cr
}


 \bold{Methods inherited from AffymetrixFile}:\cr
compareChecksum, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, readChecksum, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
