%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixLiWongModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixLiWongModel}
\docType{class}
\alias{AffymetrixLiWongModel}


 \title{The AffymetrixLiWongModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixLiWongModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixUnitGroupsModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{AffymetrixLiWongModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AffymetrixLiWongModel}\cr
extends \emph{\link[aroma.affymetrix]{AffymetrixUnitGroupsModel}}\cr



  This class represents the Li \& Wong (2001) model.
  It can be used to fit the model on a \code{\link{AffymetrixCelSet}}.
 }

 \usage{AffymetrixLiWongModel(dataSet=NULL, path=filePath("modelLiWong", getChipType(getCdf(dataSet))), model=c("pm"), ...)}

 \arguments{
   \item{...}{Arguments passed to \code{\link{AffymetrixUnitGroupsModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:fit.AffymetrixLiWongModel]{fit}} \tab Fits a unit groups (probeset) model.\cr
 \tab \code{getProbeAffinities} \tab  -\cr
}


\bold{Methods inherited from AffymetrixUnitGroupsModel}:\cr
as.character, fit, getCdf, getDataSet, getParameterSet, getPath, setup

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \section{Model estimates}{
   The estimated probe affinities are represented by the
   \code{\link{LiWongProbeAffinityFile}} class.
   Use \code{\link[aroma.affymetrix:getProbeAffinities.AffymetrixLiWongModel]{*getProbeAffinities}()} to access these.
 }

 \references{
   Li, C. and Wong, W.H. (2001), Genome Biology 2, 1-11.\cr
   Li, C. and Wong, W.H. (2001), Proc. Natl. Acad. Sci USA 98, 31-36.\cr
 }
\keyword{classes}
