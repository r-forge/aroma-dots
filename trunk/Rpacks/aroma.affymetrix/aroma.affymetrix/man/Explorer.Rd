%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Explorer.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{Explorer}
\docType{class}
\alias{Explorer}


 \title{The Explorer class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class Explorer}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{Explorer}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{ArrayExplorer}}, \link[aroma.affymetrix]{ChromosomeExplorer}\cr

public abstract static class \bold{Explorer}\cr
extends \link[R.oo]{Object}\cr


 }

 \usage{Explorer(tags="*", ...)}

 \arguments{
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be added to the output path.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{addIncludes} \tab  -\cr
 \tab \code{display} \tab  -\cr
 \tab \code{getAlias} \tab  -\cr
 \tab \code{getArrays} \tab  -\cr
 \tab \code{getAsteriskTags} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getReportPathPattern} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{getTagsOfInput} \tab  -\cr
 \tab \code{nbrOfArrays} \tab  -\cr
 \tab \code{process} \tab  -\cr
 \tab \code{setAlias} \tab  -\cr
 \tab \code{setArrays} \tab  -\cr
 \tab \code{setReportPathPattern} \tab  -\cr
 \tab \code{setSubname} \tab  -\cr
 \tab \code{setup} \tab  -\cr
 \tab \code{splitByReportPathPattern} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{Output directory structure}{
   The \emph{main directory} of an Explorer report is
    \code{reports/<name>/<subname>/}.
   The \code{<name>} is typically the same as the name of the input
   data set, and the \code{<subname>} is typically the tags of ditto.
   This main directory is where main HTML document is stored.

   For each chip type, real or "virtual" (combined), there is a
   subdirectory with the same name as the chip type, i.e.
    \code{reports/<name>/<subname>/<chiptype>/}.

   For each chip type directory, there are set of subdirectories each
   specifying a so called \emph{image layer}, e.g. an image layer
   showing the raw data, another containing the estimates of a model
   fit and so on.  Path format:
    \code{reports/<name>/<subname>/<chiptype>/<image layer>/}.
   In this directory all image files are stored, e.g. PNG files.

   In some cases one do not want to all input tags to become part of the
   subname, but instead for instance use those to name the image layer(s).
   In such cases one has to override the default names.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
 }
\keyword{classes}
