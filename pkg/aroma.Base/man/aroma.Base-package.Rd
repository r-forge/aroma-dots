%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  999.package.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{aroma.Base-package}
\alias{aroma.Base-package}
\alias{aroma.Base}
\docType{package}

\title{Package aroma.Base}


 \encoding{latin1}

 \description{
   This package provides R methods to interact with BioArray Software Environment (BASE).  The BaseFile class provides an API for easy to loading, retrieving, modifying, and saving BASE files.  The static class BasePluginDispatcher is used to call R plugins by taking care of importing data from BASE, automatic exception handling, logging to file, returning data to BASE etc.. For a reference on BASE, see [1].
   The homepage for BASE is \url{http://base.thep.lu.se/}.
 }

 \section{Requirements}{
   This package requires that the \pkg{R.oo} [2] and the \pkg{R.utils} packages are installed.  For RSP generated reports, the \pkg{R.rsp} package is required.
 }

 \section{Installation}{
   To install this package please see \url{http://www.braju.com/R/}.
 }

 \section{To get started}{
   To get started, see:
   \enumerate{
     \item \code{\link{BasePluginDispatcher}} - Generic methods to dispatcher
        R plugins for BASE. For plugin writers.  The help on
        \code{\link{main.BasePluginDispatcher}}() is highly recommended and
        explains how aroma.Base optimize memory usage etc.
     \item \code{\link{BaseFile}} - Class representing a BASE file structure.
        An instance of this class contains different types of
        \code{\link{BaseFileSection}} objects.
     \item See also the code of example plugins that come with this
        package.  They are located in
        \code{system.file("plugins", package="aroma.Base")}.
     \item For low-level reading and writing of BASE files, see
        \code{\link{readBaseFile}}() and \code{\link{writeBaseFile}}().  This is only
        for users how want to write BASE plugins without using the
        BASE plugin dispatcher, which is not recommended.
   }
 }

 \section{Additional help}{
  Below is the plugin description for the helloWorld plugin.
  Aditional plugins are available in \code{system.file("plugins", package="aroma.Base")}.
  Generic R plugin dispatcher for BASE in the aroma.Base package by Henrik Bengtsson (2005).

This specific "Hello world" plugin demonstrates the dispatchers functionality by printing a user-defined welcome message (sent from BASE), the list of arguments passed to the plugin, the code of the plugin, and then returns (here without return data).

HOW TO WRITE A PLUGIN USING AROMA.BASE:
To write an R plugin for BASE, just create a file onRun.R that defines the function

  onRun <- function(data, ..., pluginVersion="1.0") {
    # your code here
  }

which should accept a BaseFile object as the first argument, and the plugin in parameters send from BASE as additional arguments. The BASE parameters are passed to onRun() as named R arguments.  The data sent from BASE (via standard input) is already parsed and easily accessed via the methods of the BaseFile class.  See ?BaseFile in R for details.

In summary, you have to:

1) Create a plugin directory, say, ~/plugins/myPlugin/

2) Copy the generic 'basePluginDispatcher' Unix shell script to this directory. Change file rights to be executable, i.e. 'chmod ugo+x basePluginDispatcher'.  If aroma.Base is installed neither system-wide nor in '~/plugins/library/, update this file with path to local installation, e.g. add R\_LIBS="\$R\_LIBS:~hb/plugins/library".

3) Create a onRun.R defining onRun() in same directory.  All other *.R scripts in this directory will also be sourced.  R libraries installed locally in '<pluginPath>/../library/' (e.g. ~/plugins/library/) will be found automatically via 'library(<pkg>)'.

4) Setup plugin in BASE plugin settings.  Use executable '<pluginRoot>/myPlugin/basePluginDispatcher'. Data must be send as a "serial" BASEfile!

5) Test it!


REQUIREMENTS AND INSTALLATION OF AROMA.BASE:
The aroma.Base package requires R v2.1.0 or newer.  To be able to use this generic R BASE plugin dispatcher, the following R packages must be installed on the BASE system:

1) aroma.Base
2) R.utils
3) R.oo

Ideally these are installed by the system adminstrator for system-wide access, but it is also possible to install this locally.  If installed in, say, ~/plugins/library/ and the plugin is installed in ~/plugins/myPlugin/, the 'basePluginDispatcher' script will automatically be able to find the packages.  If installed locally elsewhere, modify the script according to step 2 in previous section.

To install the above packages, do the following from R:

source("http://www.braju.com/R/hbLite.R"); hbLite("aroma.Base")

To install packages locally, add argument 'lib', e.g. 'lib=~/plugins/library/'.


PROCESSING PARAMETERS WHEN LOADED:
By defining onParameters() taking plugin parameters as arguments (no 'data' object here, cf. onRun()), one can check plugin parameters as soon as they are loaded, that is, before any data is read.  Thus, if there are erroneous plugin parameters, an error can be reported immediately, without having to wait for data to be read.  Moreover, by returning modified plugin parameters in a named list, one can use onParameters() to updated parameters at this stage, leaving onRun() to focus on processing only.   Typically, one wants to convert plugin parameters, which are all character strings, to numeric values.


INTERNAL PLUGIN COMMANDS:
If plugin parameter 'internalCommands', which should contain tab-delimited strings, contains the command 'copySource', all source files are copied to the working directory. 


POST-PROCESSING SCRIPTS:
If a script named 'postPluginScript' exists in the plugin directory, it will be called (by 'basePluginDispatcher') after the R plugin finishes.  As long as this script is executable, it can be of any language, e.g. a perl or a bash script.  The data send from base is available in the 'stdin.txt' file.  The (optional) BASE data written by the R plugin is available in 'stdout.txt'.  Note that modifying 'stdout.txt' has (unfortunately) no effect!

WARNINGS AND ERRORS:
Warning and error messages are outputted to the log file as well as the standard error of the plugin dispatcher, which then in turn is imported to BASE as "plugin messages".  In addition, an error interrupts the execution of the plugin, whereas a warning does not.  Warnings and errors are \emph{not} outputted to the resulting BASE file, i.e. the standard output of the dispatcher.

TIPS'N'TRICKS:
On Windows, in order to see the 'plugin.log' in your web browser and not in Notepad, you might want to remove the file type 'LOG' from known file types. Open "This Computer", go to menu "Tools -> Folder Options..." and select tab "File Types". Find "LOG" in the list and "Remove" it. Press "OK".


BUG REPORTING:
Please report bugs to author Henrik Bengtsson <henrikb\@braju.com>.  Do not forget to report version of R, and version of aroma.Base and friends; start R and do library(aroma.Base) and cut'n'paste the startup messages.


 }

 \section{How to cite this package}{
   Whenever using this package, for now, please cite [1] as\cr

   \preformatted{@TECHREPORT{BengtssonH_2004,
  author = {Bengtsson, Henrik},
  title = {{aroma} - {A}n {R} {O}bject-oriented {M}icroarray {A}nalysis environment},
  institution = {Mathematical Statistics, Centre for Mathematical Sciences, Lund University, Sweden},
  year = {2004},
  type = {{Preprint in Mathematical Sciences}},
  number = {2004:18},
}}
\emph{}

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

 [1] L.H. Saal, C. Troein, J. Vallon-Christersson, S. Gruvberger,
     Å. Borg & C. Peterson, \emph{BioArray Software Environment (BASE):
     a platform for comprehensive management and analysis of microarray
     data}. Genome Biology, 2002.\cr

 [2] H. Bengtsson, \emph{The R.oo package - Object-Oriented Programming with References Using Standard R Code}, In Kurt Hornik, Friedrich Leisch and Achim Zeileis, editors, Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003), March 20-22, Vienna, Austria. \url{http://www.ci.tuwien.ac.at/Conferences/DSC-2003/Proceedings/}
 \cr

 [3] H. Bengtsson, \emph{A generic R plugin dispatcher for BASE},
     Poster presented at the MGED8 conference, Bergen, Norway,
     September 2005.\cr

 [4] H. Bengtsson, \emph{aroma - An R Object-oriented Microarray
     Analysis environment}, Preprints in Mathematical Sciences (manuscript
     in preparation), Mathematical Statistics, Centre for Mathematical
     Sciences, Lund University, 2004.\cr

 [5] H. Bengtsson, J. Vallon-Christersson and G. Jönsson, \emph{Calibration and assessment of channel-specific biases in microarray data with extended dynamical range}, BMC Bioinformatics, 5:177, 2004.
 \cr

 [6] Henrik Bengtsson and Ola Hössjer, \emph{Methodological Study of Affine Transformations of Gene Expression Data}, Methodological study of affine transformations of gene expression data with proposed robust non-parametric multi-dimensional normalization method, BMC Bioinformatics, 2006, 7:100.
 \cr

 }
\keyword{package}
