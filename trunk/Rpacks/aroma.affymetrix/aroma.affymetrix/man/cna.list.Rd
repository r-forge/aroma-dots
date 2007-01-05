%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  cna.R
% 
% on Sat Jan 06 01:49:26 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{cna.list}
\alias{cna.list}
\alias{cna.list}


 \title{Copy-number analysis}

 \description{
  This high-level function takes a list of probe-level data set and
  performs a complete copy-number analysis.
 }

 \usage{\method{cna}{list}(dataSets, pres="QN", plm=c("RMA", "MBEI"), posts="FLN", cnm="glad", ..., ram=1, verbose=-2)}

 \arguments{
   \item{dataSets}{A \code{\link[base]{list}} of \code{\link{AffymetrixCelSet}}.}
   \item{pres}{A \code{\link[base]{character}} \code{\link[base]{vector}} describing what methods should be used
     to preprocess the probe-level data.
     If \code{\link[base]{NULL}}, this step is skipped.}
   \item{plm}{A \code{\link[base]{character}} string describing what probe-level model (PLM)
     to fit to the preprocessed data.}
   \item{posts}{A \code{\link[base]{character}} \code{\link[base]{vector}} describing what methods should be used
     to postprocess the chip-effect estimates of the PLM.
     If \code{\link[base]{NULL}}, this step is skipped.}
   \item{cnm}{A \code{\link[base]{character}} string describing what copy-number model (CNM)
     to fit to identify copy-number regions.}
   \item{...}{Not used.}
   \item{ram}{A \code{\link[base]{double}} in (0,Inf) specifying if more or less (RAM)
     memory than default should be used.  The smaller the value is,
     the more memory conservative the analysis will be.}
   \item{verbose}{Specifies how much details should be outputted while
     analysing data.  The lower threshold, the more details.}
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{methods}
