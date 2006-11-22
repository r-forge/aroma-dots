%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenomeAnalyzer.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{GenomeAnalyzer}
\docType{class}
\alias{GenomeAnalyzer}


 \title{The GenomeAnalyzer class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class GenomeAnalyzer}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{GenomeAnalyzer}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.affymetrix]{CnPlmAnalyzer}}, \emph{\link[aroma.affymetrix]{GenomeCnAnalyzer}}, \link[aroma.affymetrix]{MbeiCnAnalyzer}, \link[aroma.affymetrix]{RmaCnAnalyzer}\cr

public abstract static class \bold{GenomeAnalyzer}\cr
extends \link[R.oo]{Object}\cr



  This abstract class provides basic methods for any type of genome analysis.
  Specific analysis, such as copy-number analysis or exon analysis, is
  done by subclasses of this class.
 }

 \usage{GenomeAnalyzer(dataSet=NULL, ...)}

 \arguments{
   \item{dataSet}{A \code{\link{AffymetrixCelSet}}.}
   \item{...}{Arguments passed to \code{\link{SnpProbeAffinityFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:fit.GenomeAnalyzer]{fit}} \tab Fits the model used by this analyzer.\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getChromosomes.GenomeAnalyzer]{getChromosomes}} \tab Gets the chromosomes available.\cr
 \tab \code{getDataSet} \tab  -\cr
 \tab \code{getGenomeInformation} \tab  -\cr
 \tab \code{\link[aroma.affymetrix:getModel.GenomeAnalyzer]{getModel}} \tab Gets the model used by this analyzer.\cr
 \tab \code{\link[aroma.affymetrix:getProgress.GenomeAnalyzer]{getProgress}} \tab Gets the fraction of units process per chromosome.\cr
 \tab \code{nbrOfArrays} \tab  -\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
