###########################################################################/**
# @RdocClass BasePluginDispatcher
#
# @title "The BasePluginDispatcher class"
#
# \description{
#  @classhierarchy
#
#  See static method \code{main()} for help.
# }
# 
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
# \section{Definitions}{
#   The \emph{plugin directory (path)} is the directory where the plugin 
#   source code and the \code{runBasePluginDispatcher} script is located. 
#   Each plugin should have a unique plugin path.
#
#   The \emph{working directory (path)} is the directory in which 
#   the plugin is evaluated by BASE.  Contrary to the plugin directory,
#   it is not possible to know the path to the working directory 
#   on before hand.
# }
#
# \section{Run the plugin dispatcher without BASE (standalone)}{
#  While developing a plugin for BASE, it is possible to test the plugin 
#  by running it on a local machine without BASE.
#  To do this, prepare the following:
#  \itemize{
#   \item Create plugin directory
#   \enumerate{
#    \item Create a directory to contain all your plugin directories, 
#          for example \code{~/plugins/}, if you have not done this already, 
#    \item Create a plugin directory, e.g. \code{~/plugins/helloWorld/}.
#    \item Copy \code{runBasePluginDispatcher} to the plugin directory.
#    \item Create a \code{onRun.R} file that defines the required \code{onRun()}
#          function.
#   }
#   \item Create a "dummy" working directory
#   \enumerate{
#    \item Create a directory to contain all your working directories, 
#          for example \code{~/work/}, if you have not done this already, 
#    \item Create a working directory, e.g. \code{~/work/12345/}.
#    \item The data sent from BASE to the plugin is available in file
#          \code{stdin.txt} in the working directory.  
#          (You can "keep" such data by selecting "Leave stdout" in the 
#          plugin settings page.)
#   }
#  }
#
#  Next is to run the plugin. Start \R in the working directory you created 
#  above, alternatively, use \code{setwd("~/work/12345/")} to change it in \R.
#  Then, in \R, do
#  \preformatted{
#   library(aroma.Base)
#   BasePluginDispatcher$main(pluginPath="~/plugins/helloWorld/")
#  }
#  This will load all \code{*.R} files in the \emph{plugin path}, read BASE 
#  data from file \code{.inputData.base}, call \code{onParameters()} as soon
#  as the plugin parameters section is read, update parameters accordingly, 
#  continue to read the rest of the data file and wrap the data up in a 
#  \code{BaseFile} object, then call then \code{onRun()} function. 
#  Any output written by the plugin (to standard output of standard error)
#  will be outputted to standard error, that is, in the example, the output
#  will be seen on the screen, cf. when running on a BASE server such output
#  will end up under "Plug-in output".
#  if valid BASE data object are returned by \code{onRun()}, they are 
#  automatically written to a BASE file named \code{.outputData.base}.
# }
#
# @author
#*/###########################################################################
setConstructorS3("BasePluginDispatcher", function(...) {
  extend(Object(), "BasePluginDispatcher")
})




#########################################################################/**
# @RdocMethod patchCode
#
# @title "Sources all R source files in the patch directory"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{pluginPath}{A @character string of the path to the plugin where 
#     the R source files for the plugin is.  For conveniency, it may also 
#     be the pathname of a file in that directory.  Typically, this should
#     be same same directory as where the calling \code{runBaseFile} 
#     script is located.  If @NULL, the the command-line argument
#     \code{--pluginPath} is used. If still @NULL, the current directory
#     is used.}
#  \item{log}{}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) a @vector of pathnames sourced. 
#   If none, @NULL is returned.
# }
#
# \section{Exception handling}{
#   All errors are caught, outputted to standard error, and ignored.
# }
#
# @author
#*/#########################################################################  
setMethodS3("patchCode", "BasePluginDispatcher", function(static, pluginPath=NULL, ..., log=FALSE) {
  # Send all output (including standard error) to the standard output.
  # This will make it possible to to use capture.output({...}) to capture
  # all output from this method.
  sink(stdout(), type="message");
  on.exit({
    sink(file=NULL, type="message"); 
  }, add=TRUE);

  # Output warnings immediately.
  oopt <- options(warn=1);
  on.exit({ options(oopt); }, add=TRUE);

  # If plugin path is not specified, get it from the command-line arguments.
  # If still not given, the current directory will be assumed.
  if (is.null(pluginPath)) {
    commandArgs <- commandArgs(asValue=TRUE, excludeReserved=TRUE)[-1];
    pluginPath <- commandArgs[["pluginPath"]];
  }

  # If pointing to, say, file aaa/bbb/ccc.txt, it is not a directory,
  # but the parent directory aaa/bbb/ is! if only 'aaa' is given
  # it might be a file, that is, not a directory; try parent again!
  if (!isDirectory(pluginPath) && isFile(pluginPath))
    pluginPath <- dirname(pluginPath);

  # The directory containing the patch code
  patchPath <- filePath(pluginPath, "..", "patch");

  if (!isDirectory(patchPath))
    return(invisible(NULL));

  # Get all *.R files in patch directory
  files <- list.files(pattern="[.]R$", path=patchPath, full.names=TRUE);
  if (length(files) == 0)
    return(invisible(NULL));

  cat("Sourcing patch code...\n");
  cat(" Patch directory contains ", length(files), 
                               " source files: ", patchPath, "\n", sep="");
 
  # Source each patch file. Errors are ignored! 
  for (file in files) {
    tryCatch({
      cat(" Sourcing patch file: ", file, "...\n", sep="");
      source(file);
      cat(" Sourcing patch file: ", file, "...done\n", sep="");
    }, error = function(ex) {
      cat(" Sourcing patch file: ", file, "...failed\n", sep="");
      print(ex);
    })
  }
  cat("Sourcing patch code...done\n");

  invisible(files);
}, static=TRUE, protected=TRUE)





#########################################################################/**
# @RdocMethod buildReports
#
# @title "Build reports from RSP templates"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{pluginPath}{A @character string of the path to the plugin where 
#     the R source files for the plugin is.  For conveniency, it may also 
#     be the pathname of a file in that directory.  Typically, this should
#     be same same directory as where the calling \code{runBaseFile} 
#     script is located.  If @NULL, the the command-line argument
#     \code{--pluginPath} is used. If still @NULL, the current directory
#     is used.}
#  \item{log}{}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) a @vector of pathnames sourced. 
#   If none, @NULL is returned.
# }
#
# \section{Exception handling}{
#   All errors are caught, outputted to standard error, and ignored.
# }
#
# @author
#*/#########################################################################  
setMethodS3("buildReports", "BasePluginDispatcher", function(static, pluginPath, envir=parent.frame(), log=FALSE, ...) {
  rspPath <- filePath(pluginPath, "rsp/");
  if (!isDirectory(rspPath)) {
    #Lc No RSP report directory found: ${rspPath}
    return(invisible(FALSE));
  }

  log && enter(log, "Generating reports from RSP templates");
  log && cat(log, "RSP directory: ", rspPath);
  log && cat(log, "Output directory: ", getwd());

  log && enter(log, "Copying the predefined contents of the main plugin include/reports/ directory to the output reports/directory");
  srcPath <- filePath(pluginPath, "..", "includes", "reports");
  if (isDirectory(srcPath)) {
    files <- copyDirectory(from=srcPath, to="reports", overwrite=TRUE);
    log && print(log, files);
  }
  log && exit(log);

  log && enter(log, "Copying the predefined contents of the plugin includes/reports/ directory to the output reports/ directory, e.g. js/ and css/ directories");
  srcPath <- filePath(rspPath, "includes", "reports");
  if (isDirectory(srcPath)) {
    files <- copyDirectory(from=srcPath, to="reports", overwrite=TRUE);
    log && print(log, files);
  }
  log && exit(log);

  msg <- capture.output({
    invisible(require(R.rsp)) || throw("Package not found/loaded: R.rsp");
  })
  log && cat(log, msg);

  log && enter(log, "Generating report chapters");
  srcPath <- filePath(rspPath, "reports");
  if (isDirectory(srcPath)) {
    sourceAllRsp(path=srcPath,
                       outputPath="reports", extension="html", envir=envir);
  }
  log && exit(log);

  log && enter(log, "Generating report index page");
  sourceAllRsp(path=rspPath, outputPath=".", extension="html", envir=envir);
  log && exit(log);

  log && exit(log);

  invisible(TRUE);
}, static=TRUE)



#########################################################################/**
# @RdocMethod main
#
# @title "Static method to start a BASE plugin"
# 
# \description{
#   @get "title".
#
#   Briefly, all R source files found in the plugin path are sourced and
#   \code{onRun()}, which should be defined in one of the files, is called.
#   Passed to \code{onRun()} is a @see "BaseFile" object for simple acccess
#   to data sent form BASE, plus the plugin parameters sent from BASE.
#   Data returned will automatically be save to file to be incorporated
#   automatically by BASE.
# }
#
# @synopsis
#
# \arguments{
#  \item{pluginPath}{A @character string of the path to the plugin where 
#     the R source files for the plugin is.  For conveniency, it may also 
#     be the pathname of a file in that directory.  Typically, this should
#     be same same directory as where the calling \code{runBaseFile} 
#     script is located.  If @NULL, the current directory is used.}
#  \item{logDetails}{A @numeric value.  The smaller (negative), the more
#     detailed the log output will be.  If zero, no logging will take
#     place.  See @see "R.utils::Verbose" for details.  The log is 
#     written to file 'plugin.log'.}
#  \item{parametersSection}{A @character @vector of regular expression
#     patterns matching the label of the BASE section in which plugin
#     parameters are stored.}
#  \item{...}{Additional arguments used for testing the plugin. 
#     For instance, using \code{stdin="stdin.txt"} will make the plugin 
#     to read BASE file data from file \code{"stdin.txt"} instead of
#     system's standard input.  Similarily, \code{stdout} and \code{stderr}
#     can be set to files.  These and other arguments are not listed in 
#     the function definition in order to minimize misuse/mistakes.
#     Experts may look at the source code for other arguments.}
# }
#
# \value{
#   Returns (invisibly) the result from \code{onRun()} as a @see "BaseFile"
#   if returned, otherwise @NULL.
# }
#
# \section{How an R BASE plugin is run}{
#  It is important to know is that BASE runs a plugin in it own unique 
#  directory, which is the current directory. All files referred to below
#  are created in this directory.
#
#  The plugin is run as follows:
#  \enumerate{
#   \item The working directory is set according to \code{workPath}. 
#         All files are created in this path.
#   \item Log message are directed to file 'plugin.log'.
#   \item A file progress bar named 'plugin.progress' is created. 
#         The size (in bytes) corresponds to 0-100\% progress.
#   \item All *.R files in the \code{pluginPath} are sourced.
#   \item BASE data is read from standard input into a @see "BaseFile"
#         object. As soon as the 'plugin parameters' section is read:
#   \enumerate{
#    \item If the plugin parameter \code{sourcePath} is set, the path
#          is sourced for additional *.R source files. 
#    \item The onParameters() function is called, allowing for early 
#          parameters validation and coersion.
#    \item The rest of the BASE file structure such as spots sections 
#          containing actually spot data is parsed/read.  In order to 
#          minimize memory usage, each spot data table is cached to a 
#          seperate file and only read into memory upon request.
#   }
#   \item Optional plugin parameters are extracted from this object.
#         The BASE file section containing the parameters is removed
#         and the resulting @see "BaseFile" object is called \code{data}.
#   \item The plugin calls
# 
#         \code{result <- onRun(data, <plugin parameters>)}
#
#         All plugin parameters are passed as name arguments, e.g. 
#         \code{pluginVersion="1.0"} etc.
#   \item When onRun(), and possibly onError(), returns, onFinally() is 
#         called with same arguments as onRun().
#   \item The \code{result} object is coerced to a @see "BaseFile" 
#          object.  If successful, it is save as a BASE file named 
#          'result.base', otherwise 'result.txt'.
#  } 
# }
#
# \section{Memory optimization}{
#   When reading the BASE file, the spot data tables for all array are
#   automatically cached to seperate files \emph{without} being read into
#   \R's memory.  The data is read into memory first when \code{getData()} 
#   is called on a @see "BaseFileSpots" object.
#   \emph{This makes it possible to pass virtually any number of arrays 
#   to an aroma.Base plugin.}
#
#   The aboves makes it extremely easy to apply array-by-array algorithms.
#   For example,
#   \preformatted{
#    onRun <- function(data, ...) {
#      # For each spot section in the BASE file
#      spotSections <- getSections(data, "spots");
#      for (spotSection in spotSections) {
#        # Get data table (automatically from temporary cached file)
#        X <- getData(spotSection, fields=signalFields);
#		   
#        # Process data
#        ...
		   
#        # Remove data from memory (but it still remains on file)
#        rm(X); gc();
#      } 
#    } # onRun()
#   }
#   When the plugin finishes, aroma.Base automatically cleans out all
#   internally cached data files.
# }
#
# \section{Exception handling}{
#  All errors are caught to avoid sudden interrupts, but the error
#  messages are sent to system's standard error.  Most of them are also
#  logged to the log file.
#
#  If an error occurs when running \code{onRun()}, it is caught and 
#  \code{onError()} is called with the caught exception \code{error} as
#  the first arguments, plus the arguments that was passed 
#  to \code{onRun()}.
#
#  Exception generated by \code{onError()} and \code{onFinally()} are
#  ignored but logged.
# }
#
# \section{Creating a new plugin}{
#  \enumerate{
#   \item Create a new plugin directory, e.g.
#                                    '~<user>/plugins/normalizeAffine/'.
#   \item Copy the \code{runBasePluginDispatcher} Unix shell script to this directory.
#   \item Copy all your *.R files to this directory, including 'onRun.R'
#         that defines \code{onRun()}.
#   \item Define the plugin on BASE plugin definition page. Use
#         '<path-to-user>/normalizeAffine/runBasePluginDispatcher' as
#         the 'execName'.
#   \item Make sure that all packages required are installed. You can 
#         install missing R packages in '~<user>/plugins/library/'.
#  }
#
#  The \code{onRun()} function must accept a @see "BaseFile" object
#  as the first argument and the optional plugin parameteters as named 
#  arguments.
#  Extra debug arguments may also be passed, which is why the function
#  should also accept these, which can be done by adding \code{...} to the
#  list of arguments. 
#  Moreover, the @see "R.utils::Verbose" object \code{log} and the
#  @see "R.utils::ProgressBar" object \code{progress} is available to
#  all onNNN() functions.
#  Example:
#
#  \preformatted{
#   onRun <- function(data, constraint=0.05, pluginVersion="1.0", ...) { 
#     log && enter(log, "Affine normalization of data");
#     on.exit(log && exit(log, "Affine normalization of data"));
#
#     # Plugin parameters are passed as strings from BASE
#     contstraint <- as.double(constraint);
#
#     # Gets reference variables to the different 'spots' sections
#     spots <- getSections(data, "spots")
#     for (spot in spots) {
#       # Get signals
#       fields <- c("intensity1", "intensity2");
#       X <- getData(spot, fields=fields);
#       X <- as.matrix(X);
#  
#       # Normalize
#       X <- normalizeAffine(X);
#  
#       # Update data
#       colnames(X) <- fields;
#       setDataFields(spot, values=X);
#     }
#
#     # Return the modified 'spots' section
#     spots;
#   } # onRun()
#  }
#
#  Note that \code{data} is already read from standard input and are very
#  easy to access via the methods of @see "BaseFile". 
#
#  For plugins that tranforms data, it is easiest to modify the actual
#  \code{data} object and return it, because it typically already have the
#  correct structure. 
# }
#
# \section{About the *.R source files}{
#  The *.R files in the plugin directory should only define
#  functions and possible some object, but not do any processing; the
#  R code should not be started until \code{onRun()} is called.
#  
#  Note that onRun() must be defined in one of the *.R files, otherwise an
#  exception is thrown saying so and the plugin fails to run. We recommend
#  that \code{onRun()} is defined in 'onRun.R'.
#  
#  Optionally, \code{onError()} and \code{onFinally()} can be defined
#  similarily.
# }
#
# \section{Install R packages locally at the plugin directory}{
#  When started, this method checks for optional \R package libraries in
#  directories named 'library/' in the \emph{parent} directory of the 
#  \code{pluginPath}.  If found, it is added at the beginning of the
#  list of package libraries that \R knowns of.  That is, one can
#  install missing/out-of-date packages in \code{<pluginPath>/../library/} 
#  so that they can be loaded in \code{onRun()}.
#
#  You can add additional library paths in one of your *.R files, e.g. 
#  in '000.R'. See \code{?.libPaths} for more details.
# }
#
# \section{About the file progress bar}{
#  Reporting the plugin progress to file by file size has the advantange
#  that it is possible to view the progress via the file system, say, by
#  ftp or ssh, and by just listing the files. No files have to be opened.
#  It is even be possible to report the progress back to the BASE user
#  interface, iff BASE would support/look for it.
# 
#  When the plugin is started, the progress is set to 0\%. When standard
#  input has been read it is 1\%, when all script files have been sourced 
#  it is 2\%, then it is up to \code{onRun()} to either call 
#  \code{increase(progress)} or \code{setProgress(progess, <fraction>)}
#  to bring the progess up to 97\% before returning. 
#  See \code{?ProgressBar} for details. When \code{onRun()} is completed, 
#  the progress is set to 98\%, when \code{onFinally()} is completed it is
#  set to 99\%. When the results has been written to file (optional), 
#  it is set to 100\%.
# }
#
# @author
#*/#########################################################################  
setMethodS3("main", "BasePluginDispatcher", function(static, pluginPath=NULL, logDetails=-99, parametersSection=c("^parameters$", "settings$"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # checkSourcePath() - Checks if a pathname is (or contains) a readable
  # directory, which may contain *.R source files.
  checkSourcePath <- function(path) {
    if (is.null(path))
      return(NULL);

    log && cat(log, level=-50, "Validating source path: ", path);

    path <- filePath(path, expandLinks="any");
    path <- getAbsolutePath(path);

    log && cat(log, level=-55, "Expanded to: ", path);

    # If pointing to, say, file aaa/bbb/ccc.txt, it is not a directory,
    # but the parent directory aaa/bbb/ is! if only 'aaa' is given
    # it might be a file, that is, not a directory; try parent again!
    if (!isDirectory(path) && isFile(path))
      path <- dirname(path);

    if (!isDirectory(path)) {
      log && cat(log, level=-50, "Not a valid source path: ", path);
      path <- NULL;
    } else {
      log && cat(log, level=-50, "Identified source path: ", path);
    }

    path;
  } # checkSourcePath()



  # Source all *.R source files in plugin directory
  loadPluginSource <- function(path, ...) {
    # Workaround: sourceDirectory() calls sourceTo() calls source() that call 
    # capabilities("iconv"), which may generate (non-catchable) message
    # 'Xlib: connection to "<host>:0.0" refused by server'; we do not 
    # need "iconv" here.
    orgCapabilities <- base::capabilities;
    assign("capabilities", function(...) FALSE, 
                                      pos=which(search() == "package:base"));
  
    tryCatch({
      # Set preprocessing hooks
      oldHooks <- getHook("sourceTo/onPreprocess");
      setHook("sourceTo/onPreprocess", function(lines, ...) { 
        lines <- LComments$compile(lines=lines);
        if (log) {
          cat(log, level=-80, "Source after pre-processing:");
          code <- displayCode(code=lines, pager="none");
          log$asGString <- FALSE;
          cat(log, level=-80, code);
          log$asGString <- TRUE;
        }
        lines;
      }, action="replace")
  
      output <- capture.output({
        files <- sourceDirectory(path=path, verbose=log, envir=pluginEnv);
      });
  
      if (length(output) > 0 && any(nchar(output) > 0)) {
        log && cat(log, "Output detected while loading plugin source:");
        log && cat(log, paste(output, collapse="\n", sep="\n"));
      }
    }, finally = {
      # Reset hooks
      setHook("sourceTo/onPreprocess", oldHooks, action="replace");
  
      # Remove above workaround
      assign("capabilities", orgCapabilities, 
                                      pos=which(search() == "package:base"));
    })
  } # loadPluginSource()



  # Source all *.R source files in plugin directory
  copyPluginSource <- function(path, dest, preprocess=FALSE, ...) {
    if (log) {
      pushState(log);
      on.exit(popState(log));
    }

    # Workaround: sourceDirectory() calls sourceTo() calls source() that call 
    # capabilities("iconv"), which may generate (non-catchable) message
    # 'Xlib: connection to "<host>:0.0" refused by server'; we do not 
    # need "iconv" here.
    orgCapabilities <- base::capabilities;
    assign("capabilities", function(...) FALSE, 
                                      pos=which(search() == "package:base"));
  
    tryCatch({
      if (preprocess) {
        # Set preprocessing hooks
        oldHooks <- getHook("sourceTo/onPreprocess");
        setHook("sourceTo/onPreprocess", function(lines, ...) { 
          lines <- LComments$compile(lines=lines);
          if (log) {
            cat(log, level=-80, "Source after pre-processing:");
            cat(log, level=-80, displayCode(code=lines, pager="none"));
          }
          lines;
        }, action="replace")
      }

      output <- capture.output({
        # Get absolute pathnames to all source files, by sourcing
        # to a dummy environment.
        pathnames <- sourceDirectory(path=path, envir=new.env());
      });

      log && enter(log, "Copying ", length(pathnames), " files");
      for (pathname in pathnames) {
        log && enter(log, "Copying file: ", pathname);

        # Get relative pathnames of the source.
        rPathname <- getRelativePath(pathname, relativeTo=path);

        # Using this, get the absolute destination pathname.
        dPathname <- filePath(dest, rPathname);

        # Make sure the path of the destination exists.
        dir <- getParent(dPathname);
        if (!isDirectory(dir)) {
          log && cat(log, "Creating directories recursively: ", dir);
          mkdirs(dir);
        }

        # Copy file
        log && cat(log, "Copying file from '", pathname, "' to '", 
                                                        dPathname, "'");
        file.copy(pathname, dPathname, overwrite=TRUE);

        log && exit(log);
      } # for(pathname ..)
      log && exit(log);
    }, finally = {
      if (preprocess) {
        # Reset hooks
        setHook("sourceTo/onPreprocess", oldHooks, action="replace");
      }
  
      # Remove above workaround
      assign("capabilities", orgCapabilities, 
                                      pos=which(search() == "package:base"));
    })
  } # copyPluginSource()



  # Calls onParameters(), which checks and modifies plugin parameters. This
  # is done as soon as they are read.  All parameters returned by 
  # onParameters() (in a list structure), will be updated and be available
  # when onRun() is called.
  onParametersRead <- function(section, ...) {
    headers <- section$headers;

    if (length(headers) == 0)
      return(NULL);

    name <- headers$section;
    found <- FALSE;
    for (pattern in parametersSection) {
      if (regexpr(pattern, name) != -1) {
        found <- TRUE;
        break;
      }
    }

    if (!found)
      return(NULL);


    # A plugin parameters section found

    # 1. The plugin parameters
    pluginArgs <- headers[names(headers) != "section"];

    log && cat(log, level=-20, "Received plugin parameters:");
    log && str(log, level=-20, pluginArgs);

    # 2. Load additional *.R source files from 'sourcePath' directory 
    #    given in the plugin parameters section.
    sourcePath <- checkSourcePath(headers$sourcePath);
    if (!is.null(sourcePath)) {
      log && enter(log, "Load addition source files from '", sourcePath, "'");
      loadPluginSource(path=sourcePath);
      log && exit(log);
    }

    # 3. Get internal plugin commands
    internalCommands <- headers$internalCommands;
    if (length(internalCommands) > 0) {
      log && cat(log, "Detected internal plugin commands: ", 
                                      paste(internalCommands, collapse=", "));
    }

    if ("copySource" %in% internalCommands) {
      if (getAbsolutePath(pluginPath) == getAbsolutePath(".")) {
        log && cat(log, "Did not copy source files, because source and destination paths are identical: ", pluginPath);
      } else {
        log && cat(log, "Copying source files: ", pluginPath);
        copyPluginSource(pluginPath, dest=".");
      }
      if (!is.null(sourcePath)) {
        if (getAbsolutePath(sourcePath) == getAbsolutePath(".")) {
          log && cat(log, "Did not copy source files, because source and destination paths are identical: ", sourcePath);
        } else {
          log && cat(log, "Copying additional source files: ", sourcePath);
          copyPluginSource(sourcePath, dest="_sourcePath");
        }
      }
    }

    # 4. Call onParameters()
    params <- NULL;
    resetWarnings();
    tryCatch({
      log && header(log, "onParameters()");
      log && enter(log, level=-20, "Calling onParameters()");
      # Assure that the correct state of the Verbose object is 
      # retained afterwards.
      log && pushState(log);
      withCallingHandlers({
        params <- eval(do.call("onParameters", pluginArgs), envir=pluginEnv);
      }, condition = function(cond) {
        log && cat(log, cond);
        log && str(log, cond);
      })
      log && popState(log);
      log && warnings(log);
      log && exit(log);
    }, error = function(ex) {
      log && popState(log);
      log && print(log, ex);
      log && warnings(log);
      log && exit(log, suffix="...failed");
      signalCondition(ex);
    })


    if (length(params) == 0)
      return(NULL);

    # 5. Check for onRun()
    onRun <- get("onRun", mode="function", envir=pluginEnv);
    if (identical(onRun, onRun.default))
      throw("No onRun() function is defined.");


    # 6. Update plugin parameters
    log && cat(log, level=-20, "Plugin parameters updated, as returned by onParameters():");
    log && str(log, params);

    # Update section headers with modified parameters.
    for (kk in seq(length=length(params))) {
      key <- names(params)[kk];
      headers[[key]] <- params[[kk]];
    }
    section$headers <- headers;

    # Return updated plugin parameter section
    list(section=section);
  } # onParametersRead()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Defaults
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Standard input
  stdin <- "stdin.txt";

  # Send all output message to standard error, since this is what BASE
  # shows on the result page as "plugin output". Do the same for error
  # message send to standard error.
  stdout <- stderr();
  stderr <- stderr();

  # Send log messages to log file
  stdlog <- "plugin.log";

  # Write results to system's standard output.
  stdres <- stdout();

  workPath <- NULL;
  baseStdin <- NULL;
  startupMessages <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Settings
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Output warnings immediately.
  oopt <- options(warn=1);
  on.exit({ options(oopt); }, add=TRUE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # "main"
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve command line argument and attach them to the
  # local environment.  This means that if the very given
  # at the command line, the will override the arguments
  # passed when calling main().
  commandArgs <- commandArgs(asValue=TRUE, excludeReserved=TRUE)[-1];
  attachLocally(commandArgs);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create a list of all arguments passed to this method
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  methodArgs <- formals(sys.function());
  excl <- c("static", "...");
  methodArgs <- methodArgs[setdiff(names(methodArgs), excl)];
  for (arg in names(methodArgs)) {
    methodArgs[[arg]] <- get(arg, inherits=FALSE);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument '...':
  args <- list(...);
  for (kk in seq(length=length(args))) {
    name <- names(args)[kk];
    value <- args[[kk]];
    assign(name, value);
    methodArgs[[name]] <- value;
    rm(name,value);
  }
  rm(args);
  

  # Argument 'stdin':
  if (is.null(stdin))
    stdin <- stdin();

  if (regexpr("\\(\\)$", stdin) != -1)
    stdin <- eval(parse(text=stdin));

  if (regexpr("\\(\\)$", stdout) != -1)
    stdout <- eval(parse(text=stdout));

  if (regexpr("\\(\\)$", stderr) != -1)
    stderr <- eval(parse(text=stderr));

  if (regexpr("\\(\\)$", stdlog) != -1)
    stdlog <- eval(parse(text=stdlog));

  if (regexpr("\\(\\)$", stdres) != -1)
    stdres <- eval(parse(text=stdres));


  # Argument 'workPath':
  if (!is.null(workPath)) {
    workPath <- Arguments$getCharacter(workPath);
    workPath <- filePath(workPath, expandLinks="any");
    if (!isDirectory(workPath))
      throw("Argument 'workPath' is not a directory: ", workPath);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set working directory. 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Note, the pluginPath is relative to the working directory, if not an
  # absolute path.
  if (!is.null(workPath)) {
    opwd <- getwd();
    setwd(workPath);
    on.exit(setwd(opwd), append=TRUE);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sink stderr 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.character(stderr)) {
    stderr <- file(stderr, open="w");
    sink(stderr, type="message");
    on.exit({
      sink(file=NULL, type="message"); 
      close(stderr);
    }, add=TRUE);
  } else if ( inherits(stderr, "connection") && 
              !identical(stderr, stderr())      ) {
    sink(stderr, type="message");
    on.exit({
      sink(file=NULL, type="message")
    }, add=TRUE);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sinking standard output to 'stdout'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if ( is.character(stdout) || (inherits(stdout, "connection") && 
                                !identical(stdout, stdout()))     ) {
    sink(stdout);
    on.exit(sink(), add=TRUE)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create log file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'logDetails':
  log <- logDetails;
  if (inherits(log, "Verbose")) {
  } else if (is.numeric(log)) {
    log <- Verbose(con=stdlog, threshold=log);
  } else {
    log <- as.logical(log);
    if (log)
      log <- Verbose(con=stdlog, threshold=-1);
  }

  if (log) {
     msg <- paste("BasePluginDispatcher v", getVersion(aroma.Base), 
                  " by ", getAuthor(aroma.Base), sep="");
     msg <- c(msg, "");
     msg <- c(msg, paste("Do library(", getName(aroma.Base), ") and ?BasePluginDispatcher in R for detailed help.", sep=""));
     header(log, msg);
     newline(log);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Details about the running system
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  log && header(log, level=-50, "R information");
  rver <- getRversion();
  rstatus <- R.Version()[["status"]];
  rdate <- paste(R.Version()[c("year", "month", "day")], collapse="-");
  msg <- sprintf("R: R v%s %s (%s)", rver, rstatus, rdate);
  info <- R.Version()[c("platform", "arch", "os", "system")];
  msg <- paste(msg, paste(names(info), ": ", info, sep="", collapse="\n"), sep="\n");
  log && cat(log, level=-50, msg);
  rm(rver, rstatus, rdate);

  log && header(log, level=-50, "System information");
  info <- .Platform;
  info <- c(info, Sys.info());
  info <- c(info, list("Current system time"=as.character(Sys.time())));
  msg <- paste(names(info), ": ", info, sep="", collapse="\n");
  log && cat(log, level=-50, msg);

  log && header(log, level=-50, "System environment variables");
  info <- as.list(Sys.getenv());
  msg <- paste(names(info), ": ", info, sep="", collapse="\n");
  log && cat(log, level=-50, msg);

  log && header(log, level=-50, "Memory information");
  info <- list();
  if (exists("memory.size", mode="function"))
    info <- c(info, memory.size=memory.size());
  if (exists("memory.limit", mode="function"))
    info <- c(info, memory.limit=memory.limit());
  if (exists("memory.size", mode="function") && 
      exists("memory.limit", mode="function")) {
    info <- c(info, memory.ratio=sprintf("%.2f%%", 
                                    100*memory.size()/memory.limit()));
  }
  msg <- paste(names(info), ": ", info, sep="", collapse="\n");
  log && cat(log, level=-50, msg);

  log && header(log, level=-50, "Locale settings");
  sysLoc <- unlist(strsplit(Sys.getlocale(), split=";"));
  names <- sub("^([^=]*)=(.*)", "\\1", sysLoc);
  sysLoc <- sub("^([^=]*)=(.*)", "\\2", sysLoc);
  names(sysLoc) <- names;
  info <- sysLoc;
  info <- c(info, Sys.localeconv());
  names <- names(info);
  names[nchar(names)==0] <- "[no name]";
  names(info) <- names;
  msg <- paste(names(info), ": ", info, sep="", collapse="\n");
  log && cat(log, level=-50, msg);
  rm(sysLoc);

  log && header(log, level=-50, "R Options");
  info <- options();
  msg <- paste(paste(names(info), ": ", info, sep="", collapse="\n"), sep="\n");
  log && cat(log, level=-50, msg);
  rm(info, msg);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Log startup messages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  log && header(log, level=-50, "Startup messages");
  if (length(startupMessages) > 0) 
    log && cat(log, level=-50, startupMessages, collapse="\n");

  log && header(log, "BasePluginDispatcher$main()");



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create a file progress bar
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  log && cat(log, level=-100, "Creating progress bar file: plugin.progess");
  progress <- FileProgressBar("plugin.progress");
  reset(progress);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Environment where plugin source is loaded and the plugin is evaluated.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # '.GlobalEnv' is given number 0 in the list of frames. Each
  #      subsequent function evaluation increases the frame stack by 1 and
  #      the call, function definition and the environment for evaluation
  #      of that function are returned by 'sys.call', 'sys.function' and
  #      'sys.frame' with the appropriate index.
  # sys.frames() gives a pairlist of all the active frames.
  # sys.nframe() returns an integer, the number of the current frame
  #      as described in the first paragraph.
  pluginEnv <- sys.frames()[[sys.nframe()]];

  # Default functions
  assign("onParameters", function(...) { list() }, envir=pluginEnv);

  onRun.default <- function(data, ...) {
    throw("Function onRun(data, ...) was not defined in any of plugin source files.");
  }

  assign("onRun", onRun.default, envir=pluginEnv);

  assign("onReport", function(data, ..., pluginPath=NULL) {
     buildReports(static, pluginPath=pluginPath, ..., log=log);
  }, envir=pluginEnv);

  assign("onError", function(data, ...) {}, envir=pluginEnv);

  assign("onFinally", function(data, ...) {}, envir=pluginEnv);

  # Make 'log' and 'progress' objects available to all onNNN() functions,
  # without having to pass them as arguments.
  assign("log", log, envir=pluginEnv);
  assign("progress", progress, envir=pluginEnv);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate the plugin path already here, if given.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The reason for validating the pluginPath here is that if it is 
  # invalid we can get an immediate error message without having to
  # read the stdin file, which may take some time.

  pluginPathOrg <- pluginPath;
  pluginPath <- checkSourcePath(pluginPath);
  if (is.null(pluginPath))
    throw("Failed to identify a valid plugin path: ", pluginPathOrg);

  # If there exists a 'library/' directory in the parent directory of
  # the plugin directory, add it to the list of library paths that R
  # knows of.
  log && enter(log, level=-50, "Checking for locally installed package libraries");
  pathname <- filePath(pluginPath, "../library/", expandLinks="any");
  pathname <- getAbsolutePath(pathname);
  if (isDirectory(pathname)) {
    # Add local library *before* global
    .libPaths(c(pathname, .libPaths()));
    log && cat(log, "Added R package library: ", pathname);
  }
  rm(pathname);
  log && cat(log, level=-10, "R library paths: ", paste(.libPaths(), collapse=", "));
  log && exit(log);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If 'pluginPath' exists, load the source now.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(pluginPath))
    loadPluginSource(path=pluginPath);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call onParameters() as soon as a parameters section has been read.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hooks <- getHook("onExit(readBaseFileSection)");
  on.exit({
    setHook("onExit(readBaseFileSection)", hooks, action="replace")
  }, add=TRUE);
  setHook("onExit(readBaseFileSection)", onParametersRead, action="append");
  log && cat(log, level=-100, "Added hook to readBaseFileSection().");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read 'stdin'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!inherits(baseStdin, "BaseFile")) {
    if (is.character(stdin)) {
      log && enter(log, level=-90, "Reading '", stdin, "'");
    } else {
      log && enter(log, level=-90, "Reading stdin");
    }

    # Read BASE file and cache all 'spots' data table on file.
    baseStdin <- BaseFile$read(stdin, extractSpotsData=TRUE);

    log && exit(log);
  } else {
    log && cat(log, "Received BASE file structure by function argument. Will not read from file or standard input.");
  }

  log && cat(log, "Validating the received BASE file structure");
  if (nbrOfSections(baseStdin) == 0) {
    throw("BaseFileException: BASE file has zero sections.");
  } else if (!isSerial(baseStdin)) {
    throw("BaseFileException: BASE file is not in a serial format.");
  } else {
    log && cat(log, level=-110, "Data received: ");
    log && print(log, level=-110, getSections(baseStdin));
  }

  setProgress(progress, 0.01);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract plugin parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in seq(parametersSection)) {
    parameters <- getSection(baseStdin, parametersSection[kk], regexpr=TRUE);
    if (!is.null(parameters)) {
      removeSection(baseStdin, parameters);
      log && cat(log, level=-110, "Data after removing 'parameters' section: ");
      log && print(log, level=-110, getSections(baseStdin));
      break;
    }
  }
  
  # Name the BaseFile object 'data' from now on.
  data <- baseStdin;
  rm(baseStdin);

  parameters <- BaseFileParameters(parameters);
  hasParameters <- (!is.null(parameters));
  if (hasParameters) {
    pluginArgs <- getParameters(parameters);
  } else {
    pluginArgs <- list();
  }

  log && cat(log, level=-20, "Data (as a list) to be passed to plugin: ");
  log && print(log, level=-20, getSections(data));

  log && cat(log, level=-20, "Plugin arguments: ");
  log && str(log, level=-20, pluginArgs);
    
  if (is.null(pluginPath)) {
    throw("Failed to identify path to plugin directory.  It is preferably given as an argument to main(), but can also be given as a BASE plugin parameter in stdin.txt.");
  }

  setProgress(progress, 0.02);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Run plugin
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  result <- NULL;

  # Call onRun()
  exception <- NULL;
  resetWarnings();
  tryCatch({
    log && header(log, "onRun()");
    log && enter(log, level=-20, "Calling onRun()");
    # About the progress bar: onRun() should increase the progress
    # from 2% to 97%, that it gets 95 steps of 1% each.
    callArgs <- c(list(data=data), pluginArgs);
    # Assure that the correct state of the Verbose object is 
    # retained afterwards.
    log && pushState(log);
    # Log warnings and immediately when they occur
    withCallingHandlers({
      result <- eval(do.call("onRun", callArgs), envir=pluginEnv);
    }, condition = function(cond) {
      log && cat(log, cond);
      log && str(log, cond);
    })
    log && popState(log);
    log && warnings(log);
    log && exit(log);
  }, error = function(ex) {
    log && popState(log);
    log && print(log, ex);
    log && warnings(log);
    log && exit(log, suffix="...failed");
    exception <<- list(exception=ex);
    tryCatch({
      log && header(log, "onError()");
      log && enter(log, level=-20, "Calling onError()");
      callArgs <- c(list(error=ex, data=data), pluginArgs);
      # Log warnings and immediately when they occur
      withCallingHandlers({
        eval(do.call("onError", callArgs), envir=pluginEnv);
      }, condition = function(cond) {
        log && cat(log, cond);
        log && str(log, cond);
      })
      log && warnings(log);
      log && exit(log);
    }, error = function(ex) {
      log && print(log, ex);
      log && warnings(log);
      log && exit(log, suffix="...failed");
    })
  })

  setProgress(progress, 0.96);


  # Call onFinally(), ignore errors
  resetWarnings();
  tryCatch({
    log && header(log, "onFinally()");
    log && enter(log, level=-20, "Calling onFinally()");
    callArgs <- c(list(data=data), pluginArgs);
    # Log warnings and immediately when they occur
    withCallingHandlers({
      eval(do.call("onFinally", callArgs), envir=pluginEnv);
    }, condition = function(cond) {
      log && cat(log, cond);
      log && str(log, cond);
    })
    log && warnings(log);
    log && exit(log);
  }, error = function(ex) {
    log && print(log, ex);
    log && warnings(log);
    log && exit(log, suffix="...failed");
  })


  log && header(log, "BasePluginDispatcher$main()");

  setProgress(progress, 0.97);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Coerce 'result' to a BaseFile object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.list(result)) {
    resetWarnings();
    tryCatch({
      log && enter(log, level=-20, "Trying to convert the result object from onRun() to a BaseFile object");
      result <- BaseFile(result);
      log && warnings(log);
      log && exit(log);
    }, error = function(ex) {
      log && print(log, ex);
      log && warnings(log);
      log && exit(log, suffix="...failed");
    })
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write 'result', if a BaseFile object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(result)) {
    cat(file=stdres, "BASEfile\nsection\tdummy\nmessage\tSince the onRun() plugin function returned nothing and BASE requires a valid BASE file as output, this dummy BASE file was generated instead in order to avoid (unnecessary) error messages.\n%\n\n");
  } else {
    if (!inherits(result, "BaseFile")) {
      log && cat(log, "WARNING: Saving the onRun() result in 'result.txt', which is not in BASE file format (and may not make sense).  This because onRun() returned an unsupported data structure:");
      log && str(log, result);
    }

    if (is.character(stdres)) {
      log && enter(log, "Writing results to file: '", stdres, "'");
    } else {
      log && enter(log, "Writing results to system's standard output");
    }

    # Writes results
    resetWarnings();
    tryCatch({
      if (inherits(result, "BaseFile")) {
        write(result, file=stdres, overwrite=TRUE);
      } else {
        write(result, file=stdres);
      }
      log && warnings(log);
      log && exit(log);
    }, error = function(ex) {
      log && print(log, ex);
      log && warnings(log);
      log && exit(log, suffix="...failed");
      exception <<- ex;
    })
  }

  setProgress(progress, 0.98);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove temporary 'spots_<assay>.txt' data files created above.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filenames <- list.files(pattern="^spots_.*[.]txt");  
  log && enter(log, level=-100, "Removing temporary data files 'spots_<assay>.txt'");
  for (file in filenames) {
    tryCatch({
      log && cat(log, level=-100, file, " [", file.info(file)$size, " bytes]");
      file.remove(file);
    }, error = function(ex) {
      log && cat(log, level=-100, "Failed to remove temporary data file (ignored): ", file);
    })
  }
  log && exit(log);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build reports, if any
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  result <- NULL;

  # Call onReport()
  exception <- NULL;
  resetWarnings();
  tryCatch({
    log && header(log, "onReport()");
    log && enter(log, level=-20, "Calling onReport()");
    callArgs <- c(list(data=data), pluginArgs, pluginPath=pluginPath);
    # Assure that the correct state of the Verbose object is 
    # retained afterwards.
    log && pushState(log);
    withCallingHandlers({
      result <- eval(do.call("onReport", callArgs), envir=pluginEnv);
    }, condition = function(cond) {
      log && cat(log, cond);
      log && str(log, cond);
    })
    log && popState(log);
    log && warnings(log);
    log && exit(log);
  }, error = function(ex) {
    log && popState(log);
    log && print(log, ex);
    log && warnings(log);
    log && exit(log, suffix="...failed");
    exception <<- list(exception=ex);
    tryCatch({
      log && header(log, "onError()");
      log && enter(log, level=-20, "Calling onError()");
      callArgs <- c(list(error=ex, data=data), pluginArgs);
      withCallingHandlers({
        eval(do.call("onError", callArgs), envir=pluginEnv);
      }, condition = function(cond) {
        log && cat(log, cond);
        log && str(log, cond);
      })
      log && warnings(log);
      log && exit(log);
    }, error = function(ex) {
      log && print(log, ex);
      log && warnings(log);
      log && exit(log, suffix="...failed");
    })

  })

  setProgress(progress, 0.99);

  if (is.null(exception)) {
    log && cat(log, "Plugin was successfully completed!");
    setProgress(progress, 1.00);
    remove(progress);
  } else {
    log && cat(log, "Plugin finished, but an exception was detected:");
    log && print(log, exception);
    setProgress(progress, 0.99);
  }

  log && cat(log, level=-50, "Giving control back to BASE. Bye!");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return stdout
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  invisible(result);
}, static=TRUE)




###########################################################################
# HISTORY:
# 2006-02-04
# o BUG FIX: memory.size() and memory.limit() do only exist on Windows.  
#   Gave an error about missing function on other platforms.
# 2006-01-13
# o Added more details on system and locale to the log file. Good for 
#   troubleshooting.
# 2005-12-15
# o Updated a grammar error in one of the log message.
# 2005-12-12
# o Added help to main() on memory optimization and how cached data tables
#   work.
# 2005-09-07
# o Now all types of conditioned are written to log file immediately by
#   the withCallingHandlers() statements.
# 2005-09-06
# o Now withCallingHandlers() is used to log warnings immediately when
#   they occur.  This will make it easier to track where warnings were
#   generated.
# o BUG FIX: Now the log is not preprocessing source code output via
#   GStrings.
# o Now the 'plugin.progress' is removed if the plugin is successful.
# o Now all copyDirectory() is automatically overwriting preexisting files.
# 2005-09-04
# o Now options() are also logged upon startup.
# o Now warnings are printed when they occur.
# o Added trial version of buildReports() and onReport().
# 2005-08-01
# o Updated the displayCode() calls to use the new argument 'code'.
# 2005-07-24
# o Now argument 'pluginPath' defaults to 'NULL' in patchCode().
# 2005-07-20
# o Renamed from BasePlugin to BasePluginDispatcher.
# o Added static patchCode() method.
# o Added Rdoc comments on how to call the plugin dispatcher standalone.
# 2005-07-08
# o After reading input BASE file, it is asserted to be in serial format.
# 2005-07-06
# o tryCatch() must not catch warnings, because then they effectively 
#   become the same as error, which interrupts the code.
# 2005-06-28
# o Added logging about "System information" upon startup.
# 2005-06-27
# o Added support for internal plugin commands via the plugin parameter
#   'internalCommands'. Currently, only the command 'copySource' is
#   supported, which copies all source files to the working directory.
# o Managed to run BASE plugins that does not give any output nor any
#   warnings or anything.
# o Now a dummy BASE file is written to stdres if onRun() did not return
#   any data.
# 2005-06-26
# o Now the result object from onRun() is written to stdout, since this
#   is where BASE expects the results to be written.
# o Pre-procesing using LComments, which are VComments with default letter
#   'L' and default name on verbose object 'log'; instead of doing
#   '#V=# log' and '#V+# <message>' one do '#L+# <message>' directly.
# o Now the plugin source is evaluated in the 'pluginEnv' environment.
# o TO DO: Load scripts in plugin path, load BASE file, and when plugin 
#   parameter section is detected, load scripts in source path, if 
#   available, and then call onParameters() and update plugin parameter 
#   section accordingly, then call onRun() and so on.
# 2005-06-23
# o Renamed 'verbose' object to 'log'; makes more sense in onRun().
# o Now pre-processing source code for VComments using sourceTo() hooks.
# o Now making use of new generic sourceDirectory() in R.utils. 
# 2005-06-20
# o Added a workaround to avoid message 'Xlib: connection to "<host>:0.0" 
#   refused by server' when calling sourceTo().
# 2005-06-19
# o Added all Rdoc comments.
# o Removed all methods but main().
# 2005-06-16
# o Created.
###########################################################################
