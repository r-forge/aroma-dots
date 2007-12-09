%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenericDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GenericDataFile}
\docType{class}
\alias{GenericDataFile}


 \title{The abstract GenericDataFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GenericDataFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{GenericDataFile}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffymetrixCdfFile}, \link[aroma.affymetrix]{AffymetrixCelFile}, \link[aroma.affymetrix]{AffymetrixCsvFile}, \link[aroma.affymetrix]{AffymetrixCsvGenomeInformation}, \link[aroma.affymetrix]{AffymetrixFile}, \link[aroma.affymetrix]{AffymetrixNetAffxCsvFile}, \link[aroma.affymetrix]{AffymetrixProbeTabFile}, \link[aroma.affymetrix]{AffymetrixTabularFile}, \link[aroma.affymetrix]{AffymetrixTsvFile}, \link[aroma.affymetrix]{AromaTabularBinaryFile}, \emph{\link[aroma.affymetrix]{AromaUflFile}}, \link[aroma.affymetrix]{AromaUgpFile}, \emph{\link[aroma.affymetrix]{AromaUnitTabularBinaryFile}}, \link[aroma.affymetrix]{ChipEffectFile}, \link[aroma.affymetrix]{CnagCfhFile}, \link[aroma.affymetrix]{CnChipEffectFile}, \link[aroma.affymetrix]{CnProbeAffinityFile}, \link[aroma.affymetrix]{DChipGenomeInformation}, \link[aroma.affymetrix]{DChipSnpInformation}, \link[aroma.affymetrix]{ExonChipEffectFile}, \link[aroma.affymetrix]{ExonProbeAffinityFile}, \link[aroma.affymetrix]{FirmaFile}, \link[aroma.affymetrix]{GenericTabularFile}, \emph{\link[aroma.affymetrix]{GenomeInformation}}, \link[aroma.affymetrix]{ParameterCelFile}, \link[aroma.affymetrix]{ProbeAffinityFile}, \link[aroma.affymetrix]{QualityAssessmentFile}, \link[aroma.affymetrix]{ResidualFile}, \link[aroma.affymetrix]{SampleAnnotationFile}, \link[aroma.affymetrix]{SnpChipEffectFile}, \emph{\link[aroma.affymetrix]{SnpInformation}}, \link[aroma.affymetrix]{SnpProbeAffinityFile}, \link[aroma.affymetrix]{UflSnpInformation}, \link[aroma.affymetrix]{UgpGenomeInformation}, \link[aroma.affymetrix]{WeightsFile}\cr

public static class \bold{GenericDataFile}\cr
extends \link[R.oo]{Object}\cr



  A GenericDataFile is an object refering to a data file on a file system.
  Note that this class is abstract and can not be instanciated, but
  instead you have to use one of the subclasses or the generic
  \code{\link[aroma.affymetrix:fromFile.GenericDataFile]{*fromFile}()} method.
 }

 \usage{GenericDataFile(filename=NULL, path=NULL, mustExist=TRUE, ...)}

 \arguments{
   \item{filename}{The filename of the file.}
   \item{path}{An optional path to the file.}
   \item{mustExist}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the file does
     not exists, otherwise not.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{compareChecksum} \tab  -\cr
 \tab \code{equals} \tab  -\cr
 \tab \code{fromFile} \tab  -\cr
 \tab \code{getAlias} \tab  -\cr
 \tab \code{getChecksum} \tab  -\cr
 \tab \code{getFilename} \tab  -\cr
 \tab \code{getFileSize} \tab  -\cr
 \tab \code{getFileType} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getPathname} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{hasTag} \tab  -\cr
 \tab \code{hasTags} \tab  -\cr
 \tab \code{isFile} \tab  -\cr
 \tab \code{readChecksum} \tab  -\cr
 \tab \code{setAlias} \tab  -\cr
 \tab \code{validateChecksum} \tab  -\cr
 \tab \code{writeChecksum} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Filename convention}{
   The filename of an \code{GenericDataFile} is structured as follows:
   \itemize{
    \item{filename}{\code{"sample001,a,b,c.CEL"}
       (this follows the \R convention (but not the Unix convention)}
    \item{fullname}{\code{"sample001,a,b,c"}}
    \item{name}{\code{"sample001"}}
    \item{tags}{\code{c("a", "b", "c")}}
    \item{extension}{\code{"CEL"}}
   }
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically part of an \code{\link{GenericDataFileSet}}.
 }
\keyword{classes}
