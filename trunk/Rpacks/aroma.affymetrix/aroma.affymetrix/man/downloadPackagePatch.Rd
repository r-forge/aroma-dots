%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  downloadPackagePatch.R
% 
% on Sun Dec 16 21:02:34 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{downloadPackagePatch}
\alias{downloadPackagePatch.default}
\alias{downloadPackagePatch}

 \title{Download a package patch}

 \description{
  Download a package patch from the package online reprocitory.
  By default, the patches are applied after being downloaded.
 }

 \usage{\method{downloadPackagePatch}{default}(pkgName, version=NULL, url=NULL, apply=TRUE, verbose=FALSE, rootPath="~/.Rpatches", ...)}

 \arguments{
   \item{pkgName}{The name of the package to be patched."}
   \item{version}{A \code{\link[base]{character}} string.}
   \item{url}{The root URL from where to download the patch.}
   \item{apply}{If \code{\link[base:logical]{TRUE}}, the patches are applied immediately after
      being downloaded.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
   \item{...}{Not used.}
 }

 \value{
   Returns (invisibly) the pathnames of the files downloaded.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link{patchPackage}}()
 }


\keyword{internal}
