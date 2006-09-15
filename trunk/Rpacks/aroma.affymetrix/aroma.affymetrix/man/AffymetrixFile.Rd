%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixFile}
\docType{class}
\alias{AffymetrixFile}


 \title{The abstract AffymetrixFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{AffymetrixFile}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffymetrixCdfFile}, \link[aroma.affymetrix]{AffymetrixCelFile}, \link[aroma.affymetrix]{ChipEffectFile}, \link[aroma.affymetrix]{CnChipEffectFile}, \link[aroma.affymetrix]{CnProbeAffinityFile}, \link[aroma.affymetrix]{GdasAnnotationFile}, \link[aroma.affymetrix]{ParameterCelFile}, \link[aroma.affymetrix]{ParameterFile}, \link[aroma.affymetrix]{ProbeAffinityFile}, \link[aroma.affymetrix]{SnpChipEffectFile}, \link[aroma.affymetrix]{SnpProbeAffinityFile}\cr

public static class \bold{AffymetrixFile}\cr
extends \link[R.oo]{Object}\cr



  An AffymetrixFile object represents a single Affymetrix file,
  e.g. an Affymetrix CEL file or an Affymetrix CDF file.
  Note that this class is abstract and can not be instanciated, but
  instead you have to use one of the subclasses or the generic
  \code{\link[aroma.affymetrix:fromFile.AffymetrixFile]{*fromFile}()} method.
 }

 \usage{AffymetrixFile(filename=NULL, path=NULL, ...)}

 \arguments{
   \item{filename}{The filename of the file.}
   \item{path}{An optional path to the file.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:as.character.AffymetrixFile]{as.character}} \tab Returns a short string describing the Affymetrix file.\cr
 \tab \code{clone} \tab  -\cr
 \tab \code{copyFile} \tab  -\cr
 \tab \code{fromFile} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getFilename.AffymetrixFile]{getFilename}} \tab Gets the filename of the file.\cr
 \tab \code{getFileSize} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getFileType.AffymetrixFile]{getFileType}} \tab Gets the file type of a file.\cr
 \tab \code{getLabel} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getName.AffymetrixFile]{getName}} \tab Gets the name of the array in the file.\cr
 \tab \code{\link[aroma.affymetrix:getPath.AffymetrixFile]{getPath}} \tab Gets the path (directory) of the file.\cr
 \tab \code{\link[aroma.affymetrix:getPathname.AffymetrixFile]{getPathname}} \tab Gets the pathname of the file.\cr
 \tab \code{setLabel} \tab  -\cr
 \tab \code{stextLabel} \tab  -\cr
 \tab \code{stextLabels} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically part of an \code{\link{AffymetrixFileSet}}.
 }
\keyword{classes}
