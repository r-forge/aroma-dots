%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  999.package.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{aroma.apd-package}
\alias{aroma.apd-package}
\alias{aroma.apd}
\docType{package}

\title{Package aroma.apd}


\description{
  The (in-house) APD file format was initially developed to store Affymetrix probe-level data, e.g. normalized CEL intensites.  Chip types can be added to APD file and similar to methods in the affxparser package, this package provides methods to read APDs organized by units (probesets).  In addition, the probe elements can be arranged optimally such that the elements are guaranteed to be read in order when, for instance, data is read unit by unit.  This speed up the read substantially.  This package is supporting other packages and not intended for interactive use.

  This package should be considered to be in "alpha" version.
}

\section{Requirements}{
  This package requires the packages \pkg{R.huge} and \pkg{affxparser}.
}

\section{Installation and updates}{

  To install this package, see \url{http://www.braju.com/R/}.
}

\section{To get started}{
  To get started, see:
  \enumerate{
    \item \code{\link{readApd}}(), \code{\link{readApdUnits}}(), \code{\link{readApdRectangle}}()
       - Reads APD files.
    \item \code{\link{celToApd}}() - creates an APD file from a CEL file.
    \item \code{\link{cdfToApdMap}}() - creates an APD read map from a CDF file.
    \item \code{\link{findApdMap}}() - finds an APD read map.
    \item \code{\link{createApd}}(), \code{\link{writeApd}}() - creates APD files.
    \item \code{\link{updateApd}}(), \code{\link{updateApdUnits}}() - updates APD files.
  }
}

\section{Search paths}{
  Typically you do not have to specify the pathname of the CDF file
  when reading APD files or similar.  Instead, the chip type is
  obtained from the APD header and the corresponding APD file is
  search for in several predefined locations.  For details how to
  specify the search path, see \code{\link[affxparser]{findCdf}}.

  Some APD files uses a corresponding read map in order to read data
  faster.  The \code{\link{findApdMap}}() method is used to find the
  corresponding APD map file containing the map vector.  How to
  specify search paths for map files, see that method.
}

\section{How to cite this package}{
  Currently no information.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\section{License}{
  The releases of this package is licensed under
  LGPL version 2.1 or newer.

  The development code of the packages is under a private licence
  (where applicable) and patches sent to the author fall under the
  latter license, but will be, if incorporated, released under the
  "release" license above.
}

\section{References}{
[1] H. Bengtsson, \emph{The R.oo package - Object-Oriented Programming with References Using Standard R Code}, In Kurt Hornik, Friedrich Leisch and Achim Zeileis, editors, Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003), March 20-22, Vienna, Austria. \url{http://www.ci.tuwien.ac.at/Conferences/DSC-2003/Proceedings/}
 \cr
}
\keyword{package}
