%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenericReporter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GenericReporter}
\docType{class}
\alias{GenericReporter}


 \title{The GenericReporter class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GenericReporter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{GenericReporter}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffymetrixCelSetReporter}}, \emph{\link[aroma.affymetrix]{AffymetrixChipTypeReporter}}, \emph{\link[aroma.affymetrix]{AffymetrixFileSetReporter}}, \emph{\link[aroma.affymetrix]{ChipTypeReporter}}, \emph{\link[aroma.affymetrix]{SpatialReporter}}\cr

public abstract static class \bold{GenericReporter}\cr
extends \link[R.oo]{Object}\cr


 }

 \usage{GenericReporter(tags="*", ...)}

 \arguments{
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be added to the output path.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getAlias} \tab  -\cr
 \tab \code{getAsteriskTags} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getInputName} \tab  -\cr
 \tab \code{getInputTags} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{process} \tab  -\cr
 \tab \code{setAlias} \tab  -\cr
 \tab \code{setup} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
