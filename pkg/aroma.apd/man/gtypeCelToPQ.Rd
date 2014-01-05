%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  gtypeCelToPQ.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{gtypeCelToPQ}
\alias{gtypeCelToPQ.default}
\alias{gtypeCelToPQ}

\title{Function to immitate Affymetrix' gtype\_cel\_to\_pq software}

\description{
 Function to immitate Affymetrix' gtype\_cel\_to\_pq software.
}

\usage{
\method{gtypeCelToPQ}{default}(filename, units=NULL, ..., cdf=NULL, nbrOfQuartets=NULL, verbose=FALSE)
}

\arguments{
  \item{filename}{The name of a CEL file.}
  \item{units}{Indices of CDF units to be returned.}
  \item{...}{Arguments passed to \code{\link[affxparser]{readCelUnits}}.}
  \item{cdf}{A CDF \code{\link[base]{list}} structure, the pathname of the CDF file, or
    \code{\link[base]{NULL}}.  If \code{\link[base]{NULL}}, the CDF file corresponding to the chip type of
    the CEL file is searched for using \code{\link[affxparser]{findCdf}}.}
  \item{nbrOfQuartets}{The number of probe quartets in the returned
    \code{\link[base]{matrix}}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an NxK \code{\link[base]{matrix}} where N is the number of probesets (SNPs) and
 K=4*Q where Q is the number of probe quartets (PMA,MMA,PMB,MMB).
 The rownames corresponds to the probeset names.
}

\examples{
# Scan for CEL files
files <- list.files(pattern="[.](c|C)(e|E)(l|L)$")

# Convert each to RAW file
for (file in files) {
  rawFile <- gsub("[.][^.]*$", ".raw", file)
  file.remove(rawFile)
  cel <- gtypeCelToPQ(file, verbose=TRUE)
  write.table(file=rawFile, cel, sep="\t", quote=FALSE)
}
}

\seealso{
 \code{\link{gtypeCelToPQ}}().
 \code{\link[affxparser]{applyCdfGroups}}.
}

\author{Henrik Bengtsson}

\references{
 [1] Affymetrix, \emph{Understanding Genotyping Probe Set Structure}, 2005.
     \url{http://www.affymetrix.com/support/developer/whitepapers/genotyping_probe_set_structure.affx}\cr
}


\keyword{programming}
