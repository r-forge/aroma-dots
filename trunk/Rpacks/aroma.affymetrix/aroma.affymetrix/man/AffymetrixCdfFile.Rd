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
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fromChipType.AffymetrixCdfFile]{fromChipType}} \tab Defines an AffymetrixCdfFile object by chip type.\cr
 \tab \code{\link[aroma.affymetrix:fromFile.AffymetrixCdfFile]{fromFile}} \tab Defines an AffymetrixCdfFile object from a CDF file.\cr
 \tab \code{\link[aroma.affymetrix:getCellIndices.AffymetrixCdfFile]{getCellIndices}} \tab Gets the cell indices unit by unit.\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getDimension} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getHeader.AffymetrixCdfFile]{getHeader}} \tab Gets the header of the CDF file.\cr
 \tab \code{\link[aroma.affymetrix:getUnitNames.AffymetrixCdfFile]{getUnitNames}} \tab Gets the names of each unit.\cr
 \tab \code{\link[aroma.affymetrix:getUnitSizes.AffymetrixCdfFile]{getUnitSizes}} \tab Gets the number of groups in each unit.\cr
 \tab \code{\link[aroma.affymetrix:indexOf.AffymetrixCdfFile]{indexOf}} \tab Gets the indices of units by their names.\cr
 \tab \code{isPm} \tab  -\cr
 \tab \code{nbrOfCells} \tab  -\cr
 \tab \code{nbrOfColumns} \tab  -\cr
 \tab \code{nbrOfQcUnits} \tab  -\cr
 \tab \code{nbrOfRows} \tab  -\cr
 \tab \code{nbrOfUnits} \tab  -\cr
 \tab \code{readUnits} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:setRestructor.AffymetrixCdfFile]{setRestructor}} \tab Specifies a function through which.\cr
}


\bold{Methods inherited from AffymetrixFile}:\cr
as.character, clone, copyFile, fromFile, getFilename, getFileSize, getFileType, getLabel, getName, getPath, getPathname, setLabel, stextLabel, stextLabels

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
