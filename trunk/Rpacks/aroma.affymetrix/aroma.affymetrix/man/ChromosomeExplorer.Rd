%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ChromosomeExplorer.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ChromosomeExplorer}
\docType{class}
\alias{ChromosomeExplorer}


 \title{The ChromosomeExplorer class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ChromosomeExplorer}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{ChromosomeExplorer}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ChromosomeExplorer}\cr
extends \link[R.oo]{Object}\cr


 }

 \usage{ChromosomeExplorer(model=NULL, tags="*", ...)}

 \arguments{
   \item{model}{A \code{\link{GladModel}} object.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be added to the output path.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{addIncludes} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:display.ChromosomeExplorer]{display}} \tab Displays the explorer in the default browser.\cr
 \tab \code{\link[aroma.affymetrix:getArrays.ChromosomeExplorer]{getArrays}} \tab Gets the names of the arrays.\cr
 \tab \code{\link[aroma.affymetrix:getChromosomes.ChromosomeExplorer]{getChromosomes}} \tab Gets the chromosomes available.\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getModel.ChromosomeExplorer]{getModel}} \tab Gets the model.\cr
 \tab \code{\link[aroma.affymetrix:getName.ChromosomeExplorer]{getName}} \tab Gets the name of the explorer.\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getTags.ChromosomeExplorer]{getTags}} \tab Gets the tags of the explorer.\cr
 \tab \code{\link[aroma.affymetrix:nbrOfArrays.ChromosomeExplorer]{nbrOfArrays}} \tab Gets the total number of arrays.\cr
 \tab \code{\link[aroma.affymetrix:process.ChromosomeExplorer]{process}} \tab Generates image files, scripts and dynamic pages for the explorer.\cr
 \tab \code{setCytoband} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Generating PNG images}{
   In order to get better looking graphs, but also to be able to generate
   bitmap images on systems without direct bitmap support, which is the case
   when running R in batch mode or on Unix without X11 support, images are
   created using the \code{\link[R.utils]{png2}} device (a wrapper for
   \code{bitmap()} immitating \code{png()}).  The \code{png()} is only
   used if \code{png2()}, which requires Ghostscript, does not.
   Note, when images are created using \code{png2()}, the images does
   not appear immediately, although the function call is completed,
   so be patient.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
  \code{\link{GladModel}}.
 }
\keyword{classes}
