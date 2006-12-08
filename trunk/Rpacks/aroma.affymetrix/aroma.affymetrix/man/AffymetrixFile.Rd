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
\link[aroma.affymetrix]{AffymetrixCdfFile}, \link[aroma.affymetrix]{AffymetrixCelFile}, \link[aroma.affymetrix]{AffymetrixProbeTabFile}, \link[aroma.affymetrix]{ChipEffectFile}, \link[aroma.affymetrix]{CnChipEffectFile}, \link[aroma.affymetrix]{CnProbeAffinityFile}, \link[aroma.affymetrix]{DChipGenomeInformation}, \link[aroma.affymetrix]{DChipSnpInformation}, \link[aroma.affymetrix]{GdasAnnotationFile}, \emph{\link[aroma.affymetrix]{GenomeInformation}}, \link[aroma.affymetrix]{GenotypeCallFile}, \link[aroma.affymetrix]{ParameterCelFile}, \link[aroma.affymetrix]{ParameterFile}, \link[aroma.affymetrix]{ProbeAffinityFile}, \link[aroma.affymetrix]{SnpChipEffectFile}, \emph{\link[aroma.affymetrix]{SnpInformation}}, \link[aroma.affymetrix]{SnpProbeAffinityFile}\cr

public static class \bold{AffymetrixFile}\cr
extends \link[R.oo]{Object}\cr



  An AffymetrixFile object represents a single Affymetrix file,
  e.g. an Affymetrix CEL file or an Affymetrix CDF file.
  Note that this class is abstract and can not be instanciated, but
  instead you have to use one of the subclasses or the generic
  \code{\link[aroma.affymetrix:fromFile.AffymetrixFile]{*fromFile}()} method.
 }

 \usage{AffymetrixFile(filename=NULL, path=NULL, mustExist=TRUE, ...)}

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
 \tab \code{\link[aroma.affymetrix:as.character.AffymetrixFile]{as.character}} \tab Returns a short string describing the Affymetrix file.\cr
 \tab \code{clone} \tab  -\cr
 \tab \code{copyFile} \tab  -\cr
 \tab \code{fromFile} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getFilename.AffymetrixFile]{getFilename}} \tab Gets the filename of the file.\cr
 \tab \code{getFileSize} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getFileType.AffymetrixFile]{getFileType}} \tab Gets the file type of a file.\cr
 \tab \code{getLabel} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getName.AffymetrixFile]{getName}} \tab Gets the name of the file.\cr
 \tab \code{\link[aroma.affymetrix:getPath.AffymetrixFile]{getPath}} \tab Gets the path (directory) of the file.\cr
 \tab \code{\link[aroma.affymetrix:getPathname.AffymetrixFile]{getPathname}} \tab Gets the pathname of the file.\cr
 \tab \code{\link[aroma.affymetrix:getTags.AffymetrixFile]{getTags}} \tab Gets the tags of the file.\cr
 \tab \code{setLabel} \tab  -\cr
 \tab \code{stextLabel} \tab  -\cr
 \tab \code{stextLabels} \tab  -\cr
 \tab \code{stextSize} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Filename convention}{
   This package assumes a special filename structure allowing methods
   to infer data-set names, sample names, file types, versions etc from
   the full path name.
 
   The filename of an \code{AffymetrixFile} is structured as follows:
   \itemize{
    \item{basename}{\code{"sample001,a,b,c.CEL"}
       (this follows the \R convention (but not the Unix convention)}
    \item{name}{\code{"sample001"}}
    \item{adjunct names}{\code{c("a", "b", "c")}}
    \item{extension}{\code{"CEL"}}
   }
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   An object of this class is typically part of an \code{\link{AffymetrixFileSet}}.
 }
\keyword{classes}
