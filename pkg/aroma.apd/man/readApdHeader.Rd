%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readApdHeader.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{readApdHeader}
\alias{readApdHeader.default}
\alias{readApdHeader}

 \title{Reads the header of an Affymetrix probe data (APD) file}

 \usage{
\method{readApdHeader}{default}(filename, ..., verbose=FALSE, .checkArgs=TRUE)
}

 \description{
   Reads the header of an Affymetrix probe data (APD) file.
 }

 \arguments{
   \item{filename}{The filename of the APD file.}
   \item{...}{Not used.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
   \item{.checkArgs}{If \code{\link[base:logical]{TRUE}}, arguments are validated, otherwise not.}
 }

 \value{
   A named \code{\link[base]{list}}.
   The \code{\link[base]{numeric}} element \code{nbrOfProbes} is the number of probe values
   available in the APD file.
   The optional \code{\link[base]{character}} element \code{name} specifies the name of
   the APD vector.
   The optional \code{\link[base]{character}} element \code{chipType} specifies the
   chip type, cf. the same field in \code{\link[affxparser]{readCelHeader}}.
   The optional \code{\link[base]{character}} element \code{maptype} specifies the type of
   probe-index map for this APD file.  Its value can be used to find
   the mapping file, see \code{\link{findApdMap}}() and \code{\link{readApdMap}}().
   All other fields are optional and \code{\link[base]{character}} values.
 }

 \details{
   The file format of an APD file is identical to the file format of an
   \code{\link[R.huge]{FileVector}}.  Most elements of the APD header are stored
   in the \code{comment} \code{\link[base]{character}} string of the file vector's header.
   The APD header \code{nbrOfProbes} is identical to the length of the
   file vector, and is \emph{not} stored in the above comment string.
 }

 \author{Henrik Bengtsson}

 \examples{\dontrun{#See ?createApd for an example.}}

 \seealso{
   \code{\link{readApd}}().
 }



\keyword{file}
\keyword{IO}
