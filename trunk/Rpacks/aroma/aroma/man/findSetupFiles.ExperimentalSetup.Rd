%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ExperimentalSetup.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findSetupFiles.ExperimentalSetup}
\alias{findSetupFiles.ExperimentalSetup}
\alias{ExperimentalSetup.findSetupFiles}
\alias{findSetupFiles.ExperimentalSetup}
\alias{findSetupFiles,ExperimentalSetup-method}


\title{Locates cDNA microarray experiment setup files}

\description{
 Locates cDNA microarray experiment setup files and returns a list of
 File object, which can be used by the static read() method in this class.
 By default it searches the current directory, but it can search any number
 of directories and recursively to an optimal depth.

 Any file with the extension \code{".setup"} is matched. The files are
 not tried to be read.
}

\synopsis{findSetupFiles.ExperimentalSetup(static, paths, recursive=TRUE, ...)}
\usage{ExperimentalSetup$findSetupFiles(paths, recursive=TRUE, ...)}

\arguments{
  \item{paths}{Vector of character strings or list of File objects
    specifying which directories to be searched.}
  \item{recursive}{If \\code{\link[base:logical]{TRUE}}, all subdirectories will also be search.
    If \\code{\link[base:logical]{FALSE}}, only the specified directory will be searched.
    If an integer is given, it specifies the maximum depth of the
    directory structure to be search; \code{recursive=0} and
    \code{recursive=FALSE} give the same results.}
}

\value{
  Returns a \code{\link[base]{list}} of filenames pointing to the setup files found.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
   files <- ExperimentalSetup$findSetupFiles(system.file("misc", package="aroma"))
   print(files)
   for (file in files) {
     setup <- ExperimentalSetup$read(file)
     print(setup)
   }
}

\seealso{
  For more information see \code{\link{ExperimentalSetup}}.
}
\keyword{internal}
\keyword{methods}
