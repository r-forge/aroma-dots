%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  999.package.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{aroma.affymetrix-package}
\alias{aroma.affymetrix-package}
\alias{aroma.affymetrix}
\docType{package}

\title{Package aroma.affymetrix}


 \description{
   This package implements classes for files and sets of files for various Affymetrix file formats, e.g AffymetrixCdfFile, AffymetrixCelFile, and AffymetrixCelSet. These are designed to be memory efficient but still being fast. The idea is to work keep all data on file and only read data into memory when needed.  Clever caching mechanisms are used to minimize the overhead of data IO.  All of the above is hidden in the package API and for the developer (and the end user), the is queried as if it lives in memory.  With this design it is only the diskspace that limits what number of arrays can be analyzed, e.g. we have quantile normalized ninety 100K SNP chips using only 150Mb RAM.  Note that the purpose of this package is to provide other developers with a solid foundation for working with Affymetrix data.  It is not intended for analysis by the end user.  Such methods should be implemented in higher-level packages utilizing aroma.affymetrix.

   This package should be considered to be in an alpha or beta phase.
   You should expect the API to be changing over time.
 }

 \section{Requirements}{
   This package requires the \pkg{aroma.apd} package and
   the \pkg{affxparser} package (Bioconductor).
   In turn, these package require packages \pkg{R.oo} [1],
   \pkg{R.utils}, and \pkg{R.huge}.
 }

 \section{Installation and updates}{

   To install this package, see instructions at
   \url{http://www.braju.com/R/}.
 }

 \section{To get started}{
   To get started, see:
   \enumerate{
     \item \code{\link{AffymetrixDataset}} - Defines a set of Affymetrix
       data files.
     \item \code{\link{AffymetrixDataFile}} - Defines an Affymetrix data
       file, e.g. CEL or APD.
   }
 }

 \section{Performance}{
   There is a performance price which we have to pay for not keeping
   data in memory but on file.  However, the performance is still
   quite good, because the underlying read methods provided by
   the \pkg{affxparser} and the \pkg{aroma.apd} packages have been
   optimized for speed.

   Note that it is much faster to access files from a local drive than
   over a local network.  Thus, you might want to consider to copy
   files on the network to a temporary local directory and work from
   there.
 }

 \section{How to cite this package}{
 }

 \section{Wishlist}{
  Here is a list of features that would be useful, but which I have
  too little time to add myself. Contributions are appreciated.
  \itemize{
    \item At the moment, nothing.
  }

  If you consider to contribute, make sure it is not already
  implemented by downloading the latest "devel" version!
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

 \references{
  Some of the reference below can be found at
  \url{http://www.maths.lth.se/bioinformatics/publications/}.\cr

 [1] H. Bengtsson, \emph{The R.oo package - Object-Oriented Programming with References Using Standard R Code}, In Kurt Hornik, Friedrich Leisch and Achim Zeileis, editors, Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003), March 20-22, Vienna, Austria. \url{http://www.ci.tuwien.ac.at/Conferences/DSC-2003/Proceedings/}
 \cr
 }
\keyword{package}
