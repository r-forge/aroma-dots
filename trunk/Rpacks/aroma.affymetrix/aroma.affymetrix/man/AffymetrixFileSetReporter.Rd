%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSetReporter.R
% 
% on Tue Nov 20 12:39:31 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixFileSetReporter}
\docType{class}
\alias{AffymetrixFileSetReporter}


 \title{The AffymetrixFileSetReporter class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixFileSetReporter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{AffymetrixFileSetReporter}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{AffymetrixCelSetReporter}}, \emph{\link[aroma.affymetrix]{SpatialReporter}}\cr

public abstract static class \bold{AffymetrixFileSetReporter}\cr
extends \link[R.oo]{Object}\cr


 }

 \usage{AffymetrixFileSetReporter(set=NULL, tags="*", ..., .setClass="AffymetrixFileSet")}

 \arguments{
   \item{set}{An \code{\link{AffymetrixFileSet}} object.}
   \item{tags}{A \code{\link[base]{character}} \code{\link[base]{vector}} of tags to be added to the output path.}
   \item{...}{Not used.}
   \item{.setClass}{The name of the class of the input set.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:getAlias.AffymetrixFileSetReporter]{getAlias}} \tab Gets the alias of the output set.\cr
 \tab \code{getAsteriskTags} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getName.AffymetrixFileSetReporter]{getName}} \tab Gets the name of the explorer.\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getTags.AffymetrixFileSetReporter]{getTags}} \tab Gets the tags of the explorer.\cr
 \tab \code{\link[aroma.affymetrix:process.AffymetrixFileSetReporter]{process}} \tab Generates report.\cr
 \tab \code{\link[aroma.affymetrix:setAlias.AffymetrixFileSetReporter]{setAlias}} \tab Sets the alias of the output set.\cr
 \tab \code{setup} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
 }
\keyword{classes}
