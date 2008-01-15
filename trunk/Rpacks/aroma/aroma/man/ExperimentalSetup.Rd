%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ExperimentalSetup.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{ExperimentalSetup}
\docType{class}
\alias{ExperimentalSetup}


\title{The ExperimentalSetup class}

\description{
 Package:  aroma \cr
\bold{Class ExperimentalSetup}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{ExperimentalSetup}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{ExperimentalSetup}\cr
extends \link[R.oo]{Object}\cr



 Represents the experimental setup of cDNA microarray experiments.
 An ExperimentalSetup object contains information about such things as
 samples used, dyes used for each sample on each hybridization, printing
 information, scanning information, which the image and data files are
 etc.
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{as.data.frame} \tab  -\cr
 \tab \code{\link[aroma:findSetupFiles.ExperimentalSetup]{findSetupFiles}} \tab Locates cDNA microarray experiment setup files.\cr
 \tab \code{getDataFiles} \tab  -\cr
 \tab \code{getDyesUsed} \tab  -\cr
 \tab \code{getImageSoftwaresUsed} \tab  -\cr
 \tab \code{getLayoutFiles} \tab  -\cr
 \tab \code{getSampleUsed} \tab  -\cr
 \tab \code{getScannersUsed} \tab  -\cr
 \tab \code{nbrOfHybridizations} \tab  -\cr
 \tab \code{print} \tab  -\cr
 \tab \code{\link[aroma:read.ExperimentalSetup]{read}} \tab Reads a file specifying the setup of a microarray experiment.\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

}


\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  \dontrun{
  # Due to a bug in R CMD check (R v1.7.1) the MicroarrayData$read() call
  # below will call getKnownSubclasses(), which will generate
  #   "Error in exists(objectName, mode = "function") :
  #	   [2003-07-07 23:32:41] Exception: F used instead of FALSE"
  # Note that the example still work, just not in R CMD check

  setup <- ExperimentalSetup$read("ScanAlyze.setup", path=system.file("misc", package="aroma"))
  print(setup)
  data <- MicroarrayData$read(getDataFiles(setup), path=system.file("misc", package="aroma"))
  raw <- getRawData(data)
  ma <- getSignal(raw)
  subplots(2*nbrOfSlides(ma))
  for (k in seq(ma)) {
    plot(ma, slide=k)
    plotSpatial(ma, slide=k)
  }
  }
}
\keyword{classes}
