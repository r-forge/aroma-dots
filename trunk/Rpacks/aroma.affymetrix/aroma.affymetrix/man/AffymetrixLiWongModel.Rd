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
\code{~~~~~~~+--}\code{\link[aroma.affymetrix]{ProbeLevelModel}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{AffymetrixLiWongModel}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.affymetrix]{AffymetrixTotalCnLiWongModel}\cr

public abstract static class \bold{AffymetrixLiWongModel}\cr
extends \emph{\link[aroma.affymetrix]{ProbeLevelModel}}\cr



  This class represents the Li \& Wong (2001) model,
  see \code{\link{AffymetrixLiWongModel}}.
 }

 \usage{AffymetrixLiWongModel(..., name="modelLiWong")}

 \arguments{
   \item{path}{The \code{\link[base]{character}} string specifying the path to the directory
      to contain the parameter-estimate files.}
   \item{...}{Arguments passed to \code{\link{ProbeLevelModel}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


\bold{Methods inherited from ProbeLevelModel}:\cr
fit, getChipEffectClass, getChipEffects, getFitFunction, getProbeAffinities, getProbeAffinityClass, plotMvsPosition

\bold{Methods inherited from AffymetrixUnitGroupsModel}:\cr
as.character, fit, getCdf, getDataSet, getLabel, getName, getParameterSet, getPath, setLabel, setup

\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
   Li, C. and Wong, W.H. (2001), Genome Biology 2, 1-11.\cr
   Li, C. and Wong, W.H. (2001), Proc. Natl. Acad. Sci USA 98, 31-36.\cr
 }
\keyword{classes}
