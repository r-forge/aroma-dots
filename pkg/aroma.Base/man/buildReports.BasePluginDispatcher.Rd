%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BasePluginDispatcher.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{buildReports.BasePluginDispatcher}
\alias{buildReports.BasePluginDispatcher}
\alias{BasePluginDispatcher.buildReports}
\alias{buildReports.BasePluginDispatcher}
\alias{buildReports,BasePluginDispatcher-method}


 \title{Build reports from RSP templates}

 \description{
   Build reports from RSP templates.
 }

 \synopsis{buildReports.BasePluginDispatcher(static, pluginPath, envir=parent.frame(), log=FALSE, ...)}
\usage{BasePluginDispatcher$buildReports(pluginPath, envir=parent.frame(), log=FALSE, ...)}

 \arguments{
  \item{pluginPath}{A \code{\link[base]{character}} string of the path to the plugin where
     the R source files for the plugin is.  For conveniency, it may also
     be the pathname of a file in that directory.  Typically, this should
     be same same directory as where the calling \code{runBaseFile}
     script is located.  If \code{\link[base]{NULL}}, the the command-line argument
     \code{--pluginPath} is used. If still \code{\link[base]{NULL}}, the current directory
     is used.}
  \item{log}{}
  \item{...}{Not used.}
 }

 \value{
   Returns (invisibly) a \code{\link[base]{vector}} of pathnames sourced.
   If none, \code{\link[base]{NULL}} is returned.
 }

 \section{Exception handling}{
   All errors are caught, outputted to standard error, and ignored.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
