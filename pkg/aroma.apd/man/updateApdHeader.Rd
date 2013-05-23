%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  updateApdHeader.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{updateApdHeader}
\alias{updateApdHeader.default}
\alias{updateApdHeader}

 \title{Updates the header of an Affymetrix probe data (APD) file}

 \usage{
\method{updateApdHeader}{default}(filename, path=NULL, ..., verbose=FALSE)
}

 \description{
   Updates the header of an Affymetrix probe data (APD) file.
 }

 \arguments{
   \item{filename}{The filename of the APD file.}
   \item{path}{The path to the APD file.}
   \item{...}{A set of named header values to be updated/added to the
      header.  A value of \code{\link[base]{NULL}} will be removed from the current header.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns (invisibly) the pathname of the file updated.
 }

 \author{Henrik Bengtsson}

 \examples{\dontrun{#See ?createApd for an example.}}

 \seealso{
   \code{\link{createApd}}() and \code{\link{updateApd}}().
 }



\keyword{file}
\keyword{IO}
