%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BasePluginDispatcher.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{BasePluginDispatcher}
\docType{class}
\alias{BasePluginDispatcher}


 \title{The BasePluginDispatcher class}

 \description{
  Package:  aroma.Base \cr
\bold{Class BasePluginDispatcher}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{BasePluginDispatcher}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BasePluginDispatcher}\cr
extends \link[R.oo]{Object}\cr



  See static method \code{main()} for help.
 }

 \usage{BasePluginDispatcher(...)}

 \arguments{
  \item{...}{Not used.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.Base:buildReports.BasePluginDispatcher]{buildReports}} \tab Build reports from RSP templates.\cr
 \tab \code{\link[aroma.Base:main.BasePluginDispatcher]{main}} \tab Static method to start a BASE plugin.\cr
 \tab \code{\link[aroma.Base:patchCode.BasePluginDispatcher]{patchCode}} \tab Sources all R source files in the patch directory.\cr
}


 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

 }

 \section{Definitions}{
   The \emph{plugin directory (path)} is the directory where the plugin
   source code and the \code{runBasePluginDispatcher} script is located.
   Each plugin should have a unique plugin path.

   The \emph{working directory (path)} is the directory in which
   the plugin is evaluated by BASE.  Contrary to the plugin directory,
   it is not possible to know the path to the working directory
   on before hand.
 }

 \section{Run the plugin dispatcher without BASE (standalone)}{
  While developing a plugin for BASE, it is possible to test the plugin
  by running it on a local machine without BASE.
  To do this, prepare the following:
  \itemize{
   \item Create plugin directory
   \enumerate{
    \item Create a directory to contain all your plugin directories,
          for example \code{~/plugins/}, if you have not done this already,
    \item Create a plugin directory, e.g. \code{~/plugins/helloWorld/}.
    \item Copy \code{runBasePluginDispatcher} to the plugin directory.
    \item Create a \code{onRun.R} file that defines the required \code{onRun()}
          function.
   }
   \item Create a "dummy" working directory
   \enumerate{
    \item Create a directory to contain all your working directories,
          for example \code{~/work/}, if you have not done this already,
    \item Create a working directory, e.g. \code{~/work/12345/}.
    \item The data sent from BASE to the plugin is available in file
          \code{stdin.txt} in the working directory.
          (You can "keep" such data by selecting "Leave stdout" in the
          plugin settings page.)
   }
  }

  Next is to run the plugin. Start \R in the working directory you created
  above, alternatively, use \code{setwd("~/work/12345/")} to change it in \R.
  Then, in \R, do
  \preformatted{
   library(aroma.Base)
   BasePluginDispatcher$main(pluginPath="~/plugins/helloWorld/")
  }
  This will load all \code{*.R} files in the \emph{plugin path}, read BASE
  data from file \code{.inputData.base}, call \code{onParameters()} as soon
  as the plugin parameters section is read, update parameters accordingly,
  continue to read the rest of the data file and wrap the data up in a
  \code{BaseFile} object, then call then \code{onRun()} function.
  Any output written by the plugin (to standard output of standard error)
  will be outputted to standard error, that is, in the example, the output
  will be seen on the screen, cf. when running on a BASE server such output
  will end up under "Plug-in output".
  if valid BASE data object are returned by \code{onRun()}, they are
  automatically written to a BASE file named \code{.outputData.base}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
