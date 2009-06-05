#########################################################################/**
# @RdocPackage aroma.Base
#
# \encoding{latin1}
#
# \description{
#   @eval "getDescription(aroma.Base)". For a reference on BASE, see [1].
#   The homepage for BASE is \url{http://base.thep.lu.se/}.
# }
#
# \section{Requirements}{
#   This package requires that the \pkg{R.oo} [2] and the \pkg{R.utils} packages are installed.  For RSP generated reports, the \pkg{R.rsp} package is required.
# }
#
# \section{Installation}{
#   To install this package please see \url{http://www.braju.com/R/}.
# }
#
# \section{To get started}{
#   To get started, see:
#   \enumerate{
#     \item @see "BasePluginDispatcher" - Generic methods to dispatcher 
#        R plugins for BASE. For plugin writers.  The help on
#        @see "main.BasePluginDispatcher" is highly recommended and
#        explains how aroma.Base optimize memory usage etc.
#     \item @see "BaseFile" - Class representing a BASE file structure.
#        An instance of this class contains different types of
#        @see "BaseFileSection" objects.
#     \item See also the code of example plugins that come with this
#        package.  They are located in 
#        \code{system.file("plugins", package="aroma.Base")}.
#     \item For low-level reading and writing of BASE files, see 
#        @see "readBaseFile" and @see "writeBaseFile".  This is only
#        for users how want to write BASE plugins without using the
#        BASE plugin dispatcher, which is not recommended.
#   }
# }
#
# \section{Additional help}{
#  Below is the plugin description for the helloWorld plugin.  
#  Aditional plugins are available in \code{system.file("plugins", package="aroma.Base")}.
#  @eval "{path <- filePath('../inst/plugins/helloWorld'); file <- list.files(pattern='^plugin.*[.]base$', path=path, full.names=TRUE)[1]; text <- readBaseFile(file)$plugin$headers$descr; text <- strsplit(text, split=c('\\\\r\\\\n', '\\\\n', '\\\\r'))[[1]]; text <- gsub('^([^a-z:]*):$', '\\\\bold{\\1}\\\\cr', text); text <- gsub(' [*]([^*]*)[*] ', ' \\\\bold{\\1} ', text); text <- gsub(' _([^_]*)_ ', ' \\\\emph{\\1} ', text); text <- gsub('([^\\\\])([_$])', '\\1\\\\\\2', text); text <- gsub('\\\\!', '!', text); text <- paste(text, collapse='\n'); text; }"
#
# }
# 
# \section{How to cite this package}{
#   Whenever using this package, for now, please cite [1] as\cr
#
#   @howtocite "aroma.Base"
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
# [1] L.H. Saal, C. Troein, J. Vallon-Christersson, S. Gruvberger, 
#     Å. Borg & C. Peterson, \emph{BioArray Software Environment (BASE): 
#     a platform for comprehensive management and analysis of microarray 
#     data}. Genome Biology, 2002.\cr
#
# [2] @include "../incl/BengtssonH_2003.bib.Rdoc" \cr
#
# [3] H. Bengtsson, \emph{A generic R plugin dispatcher for BASE}, 
#     Poster presented at the MGED8 conference, Bergen, Norway, 
#     September 2005.\cr
#
# [4] H. Bengtsson, \emph{aroma - An R Object-oriented Microarray 
#     Analysis environment}, Preprints in Mathematical Sciences (manuscript
#     in preparation), Mathematical Statistics, Centre for Mathematical
#     Sciences, Lund University, 2004.\cr
#
# [5] @include "../incl/BengtssonH_etal_2004.bib.Rdoc" \cr
#
# [6] @include "../incl/BengtssonHossjer_2006.bib.Rdoc" \cr
#
# }
#*/#########################################################################

