%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ParameterSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{ParameterSet}
\docType{class}
\alias{ParameterSet}


 \title{The abstract ParameterSet class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class ParameterSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{ParameterSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ParameterSet}\cr
extends \link[R.oo]{Object}\cr



  An ParameterSet object represents a set of parameter estimates.
 }

 \usage{ParameterSet(files=NULL, ...)}

 \arguments{
   \item{filename}{The filename of the data file.}
   \item{path}{An optional path to the data data file.}
   \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.affymetrix:000-LB-.ParameterSet]{[}} \tab Extract a subset of the parameter set.\cr
 \tab \code{\link[aroma.affymetrix:000-LB--LB-.ParameterSet]{[[}} \tab Gets a subset of the files in the parameter set.\cr
 \tab \code{\link[aroma.affymetrix:as.character.ParameterSet]{as.character}} \tab Returns a short string describing this parameter set.\cr
 \tab \code{\link[aroma.affymetrix:as.list.ParameterSet]{as.list}} \tab Returns the data files of the parameter set.\cr
 \tab \code{\link[aroma.affymetrix:extract.ParameterSet]{extract}} \tab Extract a subset of the parameter set.\cr
 \tab \code{\link[aroma.affymetrix:fromFiles.ParameterSet]{fromFiles}} \tab Defines a ParameterSet object by searching for files.\cr
 \tab \code{\link[aroma.affymetrix:getName.ParameterSet]{getName}} \tab Gets the name of this parameter set.\cr
 \tab \code{\link[aroma.affymetrix:getPath.ParameterSet]{getPath}} \tab Gets the path (directory) of this parameter set.\cr
 \tab \code{\link[aroma.affymetrix:lapply.ParameterSet]{lapply}} \tab Applies a function to each of the parameter files.\cr
 \tab \code{\link[aroma.affymetrix:nbrOfArrays.ParameterSet]{nbrOfArrays}} \tab Gets the number of arrays in this parameter set.\cr
 \tab \code{\link[aroma.affymetrix:seq.ParameterSet]{seq}} \tab Gets an vector of file indices.\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson \url{http://www.braju.com/R/}}

\keyword{classes}
\keyword{internal}
