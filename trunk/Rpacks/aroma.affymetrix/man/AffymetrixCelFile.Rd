%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixCelFile}
\docType{class}
\alias{AffymetrixCelFile}


 \title{The AffymetrixCelFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixCelFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{AffymetrixCelFile}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{LiWongProbeAffinityFile}, \emph{\link[aroma.affymetrix]{ParameterCelFile}}\cr

public static class \bold{AffymetrixCelFile}\cr
extends \link[aroma.affymetrix]{AffymetrixFile}\cr



  An AffymetrixCelFile object represents a single Affymetrix CEL file.
 }

 \usage{AffymetrixCelFile(..., cdf=NULL)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixFile}}.}
   \item{cdf}{An optional \code{\link{AffymetrixCdfFile}} making it possible to
     override the default CDF file as specified by the CEL file header.
     The requirement is that its number of cells must match that of
     the CEL file.
     If \code{\link[base]{NULL}}, the CDF structure is inferred from the the chip type
     as specified in the CEL file header.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:fromFile.AffymetrixCelFile]{fromFile}} \tab Defines an AffymetrixCelFile object from a CEL file.\cr
 \tab \code{\link[aroma.affymetrix:getCdf.AffymetrixCelFile]{getCdf}} \tab Gets the CDF structure for this CEL file.\cr
 \tab \code{setCdf} \tab  -\cr
}


\bold{Methods inherited from AffymetrixFile}:\cr
as.character, fromFile, getFilename, getFileSize, getFileType, getName, getPath, getPathname

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{For developers}{
   If you subclass this class, please make sure to query the
   \code{\link{AffymetrixCdfFile}} object (see \code{\link[aroma.affymetrix:getCdf.AffymetrixCelFile]{*getCdf}()}) whenever
   querying CDF information.  Do not use the CDF file inferred from the
   chip type in CEL header, unless you really want it to be hardwired that
   way, otherwise you will break to possibility to override the CDF
   structure.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
