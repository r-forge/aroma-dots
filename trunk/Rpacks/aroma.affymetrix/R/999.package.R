#########################################################################/**
# @RdocPackage aroma.affymetrix
#
# \description{
#   @eval "getDescription(aroma.affymetrix)" 
#
#   This package should be considered to be in an alpha or beta phase.
#   You should expect the API to be changing over time.
# }
#
# \section{Requirements}{
#   This package requires the \pkg{aroma.apd} package and 
#   the \pkg{affxparser} package (Bioconductor).
#   In turn, these package require packages \pkg{R.oo} [1], 
#   \pkg{R.utils}, and \pkg{R.huge}.
# }
#
# \section{Installation and updates}{
#
#   To install this package, see instructions at 
#   \url{http://www.braju.com/R/}.
# } 
#
# \section{To get started}{
#   To get started, see:
#   \enumerate{
#     \item @see "AffymetrixDataset" - Defines a set of Affymetrix 
#       data files.
#     \item @see "AffymetrixDataFile" - Defines an Affymetrix data 
#       file, e.g. CEL or APD.
#   }
# }
# 
# \section{Performance}{
#   There is a performance price which we have to pay for not keeping 
#   data in memory but on file.  However, the performance is still
#   quite good, because the underlying read methods provided by
#   the \pkg{affxparser} and the \pkg{aroma.apd} packages have been 
#   optimized for speed.
#
#   Note that it is much faster to access files from a local drive than
#   over a local network.  Thus, you might want to consider to copy
#   files on the network to a temporary local directory and work from
#   there.
# }
#
# \section{How to cite this package}{
# }
#
# \section{Wishlist}{
#  Here is a list of features that would be useful, but which I have
#  too little time to add myself. Contributions are appreciated.
#  \itemize{
#    \item At the moment, nothing.
#  }
#
#  If you consider to contribute, make sure it is not already 
#  implemented by downloading the latest "devel" version!
# } 
#
# @author
#
# \section{License}{
#   The releases of this package is licensed under 
#   LGPL version 2.1 or newer.
#
#   The development code of the packages is under a private licence 
#   (where applicable) and patches sent to the author fall under the
#   latter license, but will be, if incorporated, released under the
#   "release" license above.
# }
# 
# \references{
#  Some of the reference below can be found at 
#  \url{http://www.maths.lth.se/bioinformatics/publications/}.\cr
#
# [1] @include "../incl/BengtssonH_2003.bib.Rdoc" \cr
# }
#*/#########################################################################

