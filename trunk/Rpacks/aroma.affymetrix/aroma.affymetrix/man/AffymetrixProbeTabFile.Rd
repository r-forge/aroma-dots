%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixProbeTabFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AffymetrixProbeTabFile}
\docType{class}
\alias{AffymetrixProbeTabFile}


 \title{The AffymetrixProbeTabFile class}

 \description{
  Package:  aroma.affymetrix \cr
\bold{Class AffymetrixProbeTabFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.affymetrix]{AffymetrixFile}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{AffymetrixProbeTabFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{AffymetrixProbeTabFile}\cr
extends \link[aroma.affymetrix]{AffymetrixFile}\cr



  An AffymetrixProbeTabFile represents an interface to query the data
  contained in an Affymetrix probe tab file, e.g.
  \code{Mapping250K_Nsp_probe_tab}.
 }

 \usage{AffymetrixProbeTabFile(...)}

 \arguments{
   \item{...}{Arguments passed to the constructor of \code{\link{AffymetrixFile}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fromCdf} \tab  -\cr
 \tab \code{fromChipType} \tab  -\cr
 \tab \code{getCdf} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
}


 \bold{Methods inherited from AffymetrixFile}:\cr
compareChecksum, fromFile, getAlias, getAttribute, getAttributes, getChecksum, getFilename, getFileSize, getFileType, getFullName, getName, getPath, getPathname, getTags, hasTag, hasTags, readChecksum, setAlias, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, testAttributes, validateChecksum, writeChecksum

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


 }

 \section{About probe-tab files}{
  Probe-tab files are provided by Affymetrix and contains information
  on the probes.  Note that not necessarily all probes are represented
  in the files, e.g. typically only PM probes are given and MM are
  left to be inferred from the PMs.

  The below is an extract of the \code{Mapping250K_Nsp_probe_tab} file
  obtained from the Affymetrix website.  Note that columns are separated
  by TABs.
  \preformatted{
   SNP_A-1780270	1633	2398	3	TTGTTAAGCAAGTGAGTTATTTTAT	f	PM	C
   SNP_A-1780270	1633	2399	3	TTGTTAAGCAAGTGACTTATTTTAT	f	PM	G
   SNP_A-1780270	1951	1780	-4	GGATAAAATAAAATAACTCACTTGC	r	PM	C
   ...
   SNP_A-4241299	2553	1658	4	AAACACATTTTTGGGTCGTAAGGAA	f	PM	G
  }
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{classes}
\keyword{internal}
