%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BasePluginDispatcher.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{main.BasePluginDispatcher}
\alias{main.BasePluginDispatcher}
\alias{BasePluginDispatcher.main}
\alias{main.BasePluginDispatcher}
\alias{main,BasePluginDispatcher-method}


 \title{Static method to start a BASE plugin}

 \description{
   Static method to start a BASE plugin.

   Briefly, all R source files found in the plugin path are sourced and
   \code{onRun()}, which should be defined in one of the files, is called.
   Passed to \code{onRun()} is a \code{\link{BaseFile}} object for simple acccess
   to data sent form BASE, plus the plugin parameters sent from BASE.
   Data returned will automatically be save to file to be incorporated
   automatically by BASE.
 }

 \synopsis{main.BasePluginDispatcher(static, pluginPath=NULL, logDetails=-99, parametersSection=c("^parameters$", "settings$"), ...)}
\usage{BasePluginDispatcher$main(pluginPath=NULL, logDetails=-99, parametersSection=c("^parameters$", "settings$"), ...)}

 \arguments{
  \item{pluginPath}{A \code{\link[base]{character}} string of the path to the plugin where
     the R source files for the plugin is.  For conveniency, it may also
     be the pathname of a file in that directory.  Typically, this should
     be same same directory as where the calling \code{runBaseFile}
     script is located.  If \code{\link[base]{NULL}}, the current directory is used.}
  \item{logDetails}{A \code{\link[base]{numeric}} value.  The smaller (negative), the more
     detailed the log output will be.  If zero, no logging will take
     place.  See \code{\link[R.utils]{Verbose}} for details.  The log is
     written to file 'plugin.log'.}
  \item{parametersSection}{A \code{\link[base]{character}} \code{\link[base]{vector}} of regular expression
     patterns matching the label of the BASE section in which plugin
     parameters are stored.}
  \item{...}{Additional arguments used for testing the plugin.
     For instance, using \code{stdin="stdin.txt"} will make the plugin
     to read BASE file data from file \code{"stdin.txt"} instead of
     system's standard input.  Similarily, \code{stdout} and \code{stderr}
     can be set to files.  These and other arguments are not listed in
     the function definition in order to minimize misuse/mistakes.
     Experts may look at the source code for other arguments.}
 }

 \value{
   Returns (invisibly) the result from \code{onRun()} as a \code{\link{BaseFile}}
   if returned, otherwise \code{\link[base]{NULL}}.
 }

 \section{How an R BASE plugin is run}{
  It is important to know is that BASE runs a plugin in it own unique
  directory, which is the current directory. All files referred to below
  are created in this directory.

  The plugin is run as follows:
  \enumerate{
   \item The working directory is set according to \code{workPath}.
         All files are created in this path.
   \item Log message are directed to file 'plugin.log'.
   \item A file progress bar named 'plugin.progress' is created.
         The size (in bytes) corresponds to 0-100\% progress.
   \item All *.R files in the \code{pluginPath} are sourced.
   \item BASE data is read from standard input into a \code{\link{BaseFile}}
         object. As soon as the 'plugin parameters' section is read:
   \enumerate{
    \item If the plugin parameter \code{sourcePath} is set, the path
          is sourced for additional *.R source files.
    \item The onParameters() function is called, allowing for early
          parameters validation and coersion.
    \item The rest of the BASE file structure such as spots sections
          containing actually spot data is parsed/read.  In order to
          minimize memory usage, each spot data table is cached to a
          seperate file and only read into memory upon request.
   }
   \item Optional plugin parameters are extracted from this object.
         The BASE file section containing the parameters is removed
         and the resulting \code{\link{BaseFile}} object is called \code{data}.
   \item The plugin calls

         \code{result <- onRun(data, <plugin parameters>)}

         All plugin parameters are passed as name arguments, e.g.
         \code{pluginVersion="1.0"} etc.
   \item When onRun(), and possibly onError(), returns, onFinally() is
         called with same arguments as onRun().
   \item The \code{result} object is coerced to a \code{\link{BaseFile}}
          object.  If successful, it is save as a BASE file named
          'result.base', otherwise 'result.txt'.
  }
 }

 \section{Memory optimization}{
   When reading the BASE file, the spot data tables for all array are
   automatically cached to seperate files \emph{without} being read into
   \R's memory.  The data is read into memory first when \code{getData()}
   is called on a \code{\link{BaseFileSpots}} object.
   \emph{This makes it possible to pass virtually any number of arrays
   to an aroma.Base plugin.}

   The aboves makes it extremely easy to apply array-by-array algorithms.
   For example,
   \preformatted{
    onRun <- function(data, ...) {
      # For each spot section in the BASE file
      spotSections <- getSections(data, "spots");
      for (spotSection in spotSections) {
        # Get data table (automatically from temporary cached file)
        X <- getData(spotSection, fields=signalFields);
		 
        # Process data
        ...
        # Remove data from memory (but it still remains on file)
        rm(X); gc();
      }
    } # onRun()
   }
   When the plugin finishes, aroma.Base automatically cleans out all
   internally cached data files.
 }

 \section{Exception handling}{
  All errors are caught to avoid sudden interrupts, but the error
  messages are sent to system's standard error.  Most of them are also
  logged to the log file.

  If an error occurs when running \code{onRun()}, it is caught and
  \code{onError()} is called with the caught exception \code{error} as
  the first arguments, plus the arguments that was passed
  to \code{onRun()}.

  Exception generated by \code{onError()} and \code{onFinally()} are
  ignored but logged.
 }

 \section{Creating a new plugin}{
  \enumerate{
   \item Create a new plugin directory, e.g.
                                    '~<user>/plugins/normalizeAffine/'.
   \item Copy the \code{runBasePluginDispatcher} Unix shell script to this directory.
   \item Copy all your *.R files to this directory, including 'onRun.R'
         that defines \code{onRun()}.
   \item Define the plugin on BASE plugin definition page. Use
         '<path-to-user>/normalizeAffine/runBasePluginDispatcher' as
         the 'execName'.
   \item Make sure that all packages required are installed. You can
         install missing R packages in '~<user>/plugins/library/'.
  }

  The \code{onRun()} function must accept a \code{\link{BaseFile}} object
  as the first argument and the optional plugin parameteters as named
  arguments.
  Extra debug arguments may also be passed, which is why the function
  should also accept these, which can be done by adding \code{...} to the
  list of arguments.
  Moreover, the \code{\link[R.utils]{Verbose}} object \code{log} and the
  \code{\link[R.utils]{ProgressBar}} object \code{progress} is available to
  all onNNN() functions.
  Example:

  \preformatted{
   onRun <- function(data, constraint=0.05, pluginVersion="1.0", ...) {
     log && enter(log, "Affine normalization of data");
     on.exit(log && exit(log, "Affine normalization of data"));

     # Plugin parameters are passed as strings from BASE
     contstraint <- as.double(constraint);

     # Gets reference variables to the different 'spots' sections
     spots <- getSections(data, "spots")
     for (spot in spots) {
       # Get signals
       fields <- c("intensity1", "intensity2");
       X <- getData(spot, fields=fields);
       X <- as.matrix(X);

       # Normalize
       X <- normalizeAffine(X);

       # Update data
       colnames(X) <- fields;
       setDataFields(spot, values=X);
     }

     # Return the modified 'spots' section
     spots;
   } # onRun()
  }

  Note that \code{data} is already read from standard input and are very
  easy to access via the methods of \code{\link{BaseFile}}.

  For plugins that tranforms data, it is easiest to modify the actual
  \code{data} object and return it, because it typically already have the
  correct structure.
 }

 \section{About the *.R source files}{
  The *.R files in the plugin directory should only define
  functions and possible some object, but not do any processing; the
  R code should not be started until \code{onRun()} is called.

  Note that onRun() must be defined in one of the *.R files, otherwise an
  exception is thrown saying so and the plugin fails to run. We recommend
  that \code{onRun()} is defined in 'onRun.R'.

  Optionally, \code{onError()} and \code{onFinally()} can be defined
  similarily.
 }

 \section{Install R packages locally at the plugin directory}{
  When started, this method checks for optional \R package libraries in
  directories named 'library/' in the \emph{parent} directory of the
  \code{pluginPath}.  If found, it is added at the beginning of the
  list of package libraries that \R knowns of.  That is, one can
  install missing/out-of-date packages in \code{<pluginPath>/../library/}
  so that they can be loaded in \code{onRun()}.

  You can add additional library paths in one of your *.R files, e.g.
  in '000.R'. See \code{?.libPaths} for more details.
 }

 \section{About the file progress bar}{
  Reporting the plugin progress to file by file size has the advantange
  that it is possible to view the progress via the file system, say, by
  ftp or ssh, and by just listing the files. No files have to be opened.
  It is even be possible to report the progress back to the BASE user
  interface, iff BASE would support/look for it.

  When the plugin is started, the progress is set to 0\%. When standard
  input has been read it is 1\%, when all script files have been sourced
  it is 2\%, then it is up to \code{onRun()} to either call
  \code{increase(progress)} or \code{setProgress(progess, <fraction>)}
  to bring the progess up to 97\% before returning.
  See \code{?ProgressBar} for details. When \code{onRun()} is completed,
  the progress is set to 98\%, when \code{onFinally()} is completed it is
  set to 99\%. When the results has been written to file (optional),
  it is set to 100\%.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
