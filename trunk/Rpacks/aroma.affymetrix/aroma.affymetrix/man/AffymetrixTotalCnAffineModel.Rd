%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixTotalCnAffineModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixTotalCnAffineModel}
\docType{class}
\alias{AffymetrixTotalCnAffineModel}


 \title{The AffymetrixTotalCnAffineModel class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixTotalCnAffineModel}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixUnitGroupsModel}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.affymetrix]{AffymetrixAffineModel}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{AffymetrixTotalCnAffineModel}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AffymetrixTotalCnAffineModel}\cr
extends \emph{\link[aroma.affymetrix]{AffymetrixAffineModel}}\cr



  This class represents the affine model used by Bengtsson \& Hössjer (2006).
  It can be used to fit the model on a \code{\link{AffymetrixCelSet}}.
 }

 \usage{AffymetrixTotalCnAffineModel(..., name="modelTotalCnAffine")}

 \arguments{
   \item{dataSet}{An \code{\link{AffymetrixCelSet}} object.}
   \item{path}{The \code{\link[base]{character}} string specifying the path to the directory
      to contain the parameter-estimate files.}
   \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


\bold{Methods inherited from AffymetrixAffineModel}:\cr
getFitFunction, getProbeAffinityClass

\bold{Methods inherited from ProbeLevelModel}:\cr
fit, getChipEffectClass, getChipEffects, getFitFunction, getProbeAffinities, getProbeAffinityClass, highlight, plotMvsPosition

\bold{Methods inherited from AffymetrixUnitGroupsModel}:\cr
as.character, fit, getCdf, getDataSet, getLabel, getName, getParameterSet, getPath, setLabel, setup

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \section{Model estimates}{
   The estimated probe affinities are represented by the
   \code{\link{AffineProbeAffinityFile}} class.
 }

 \references{
 }
\keyword{classes}
