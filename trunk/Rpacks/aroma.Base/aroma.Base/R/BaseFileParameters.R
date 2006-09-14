###########################################################################/**
# @RdocClass BaseFileParameters
#
# @title "The BaseFileParameters class extending BaseFileSection"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{section}{A section @list structure.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
# \details{
# }
#
#
# @author
#
# \seealso{
#   @see "BaseFile".
# }
#*/###########################################################################
setConstructorS3("BaseFileParameters", function(section=NULL, ...) {
  if (inherits(section, "BaseFileParameters"))
    return(clone(section));

  if (inherits(section, "BaseFileSection")) {
    section <- as.list(section);
    section$headers$section <- "parameters";
  }

  if (is.list(section)) {
    headers <- section$headers;
    if (length(headers) == 0)
      throw("Argument 'section' has not or empty 'headers' element.");

    if (headers$section != "parameters")
      throw("Not a 'parameters' section.");
  } else if (!is.null(section)) {
    throw("Argument 'section' is not a list nor a BaseFileSection object: ", 
                                                             class(section));
  }

  extend(BaseFileSection(section, type="parameters"), "BaseFileParameters")
})




#########################################################################/**
# @RdocMethod as.character
#
# @title "Gets a string description of BASE parameters"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################  
setMethodS3("as.character", "BaseFileParameters", function(this, ...) {
  s <- paste(class(this)[1], ": BASE section '", this$type, "'.", sep="");
  if (nbrOfParameters(this) > 0) {
    s <- paste(s, " ", nbrOfParameters(this), " parameters: ", paste(names(getParameters(this)), collapse=", "), ".", sep="");
  }
  s;
})



#########################################################################/**
# @RdocMethod nbrOfParameters
#
# @title "Gets the number of parameters in a BASE parameters section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################  
setMethodS3("nbrOfParameters", "BaseFileParameters", function(this, ...) {
  length(getHeaders(this));
})



#########################################################################/**
# @RdocMethod attachParameters
#
# @title "Assigns the parameters of a BASE parameters section locally"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getParameters".}
#  \item{envir}{The @environment where to assign the parameters.}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################  
setMethodS3("attachParameters", "BaseFileParameters", function(this, ..., envir=parent.frame()) {
  parameters <- getParameters(this, ...);
  for (parameter in parameters) {
    assign(parameter, value, envir=envir);
  }
})



#########################################################################/**
# @RdocMethod getParameters
#
# @title "Gets parameters from a BASE parameters section"
# 
# \description{
#   @get "title".  All parameters are retrieved via 
#   @seemethod "getParameter" so that certain parameters are pre-processed.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getHeader".}
# }
#
# \value{
#   Returns a named @list structure of parameters.
# }
#
# @author
#
# \seealso{
#   @seemethod "getParameter".
#   @seeclass
# }
#*/#########################################################################  
setMethodS3("getParameters", "BaseFileParameters", function(this, ...) {
  headers <- getHeaders(this, ...);
  # Pre-process all parameters via getParameter()
  res <- list();
  for (name in names(headers))
    res[[name]] <- getParameter(this, name);
  res;
})



#########################################################################/**
# @RdocMethod getParameter
#
# @title "Gets one parameter from a BASE parameters section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{name}{A @character string of the name of the parameter.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns the value of a parameter.
# }
#
# \details{
#   Certain parameters are pre-processed. This is a trial behavior. See
#   source code, i.e. \code{print(getParameter.BaseFileParameters)} 
#   for details.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################  
setMethodS3("getParameter", "BaseFileParameters", function(this, name, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  safeParse <- function(code, ...) {
    if (is.null(code))
      return(NULL);

    tryCatch({
      expr <- parse(text=code);
    }, error = function(ex) {
      throw("Invalid code '", code, "'. Reason was: ", getMessage(ex));
    })

    expr;
  } # safeParse()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # "main"
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  value <- getHeader(this, name, ...);

  # Pre-process certain parameters
  if (name == "language") {
    # Default value is "R"
    if (is.null(value))
      value <- "R";
  } else if (name == "pluginPath") {
    value <- filePath(value, expandLinks="any");
  } else if (name == "setupCode") {
    value <- safeParse(value);
  } else if (name == "runCode") {
    value <- safeParse(value);
  } else if (name == "onRun") {
    if (!is.null(value)) {
      if (!exists(value, mode="function"))
        throw("Invalid onRun parameter; no such function: ", value);
      value <- exists(value, mode="function");
    }
  }

  value;
})



#########################################################################/**
# @RdocMethod getPluginVersion
#
# @title "Gets the plugin version of a BASE parameters section"
# 
# \description{
#   @get "title".  
#   This parameter is automatically passed by BASE for transformation
#   plugins.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string, or @NULL, if missing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################  
setMethodS3("getPluginVersion", "BaseFileParameters", function(this, ...) {
  value <- getHeader(this, "pluginVersion", ...);
  if (is.null(value))
    return(NULL);
  as.character(value);
})



############################################################################
# HISTORY: 
# 2005-06-19
# o Added all Rdoc comments.
# 2005-06-16
# o Created.
############################################################################
