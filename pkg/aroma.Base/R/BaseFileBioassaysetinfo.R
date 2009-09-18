###########################################################################/**
# @RdocClass BaseFileBioassaysetinfo
#
# @title "The BaseFileBioassaysetinfo class extending BaseFileSection"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to constructor of superclass.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
# @author
#*/###########################################################################
setConstructorS3("BaseFileBioassaysetinfo", function(...) {
  extend(BaseFileSection(..., type="bioassaysetinfo"), "BaseFileBioassaysetinfo")
})



#########################################################################/**
# @RdocMethod as.character
#
# @title "Gets a string description of object"
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
setMethodS3("as.character", "BaseFileBioassaysetinfo", function(x, ...) {
  # To please R CMD check
  this <- x;
  getLayout <- NULL; rm(getLayout);

  s <- paste(class(this)[1], ":", sep="");
  s <- paste(s, " Type: ", getType(this), sep="");
  assays <- getName(this);
  s <- paste(s, ". Number of assays: ", nrow(assays), sep="");
  if (length(assays) <= 4) {
    assays <- paste(assays, collapse=", ");
  } else {
    assays <- c(assays[c(1,2)], "...", assays[length(assays)]);
  }
  s <- paste(s, " (", paste(assays, collapse=", "), ")", sep="");
  tryCatch({
    layout <- getLayout(this);
    s <- paste(s, ". ", as.character(layout), sep="");
  }, error = function(ex) {});
  s <- paste(s, ". Number of headers: ", nbrOfHeaders(this), sep="");
  s;
})



#########################################################################/**
# @RdocMethod validate
#
# @title "Validates bioassay-set information section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{action}{A @character string specifying action if a problem is 
#               detected. See @see "validate.BaseFileSection" for details.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string explaining problems detected. If no problems
#   where detected, @NULL is returned.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("validate", "BaseFileBioassaysetinfo", function(this, action=c("error", "warning", "quiet"), ...) {
  action <- match.arg(action);

  msg <- NextMethod("validate", action="none");

  name <- getName(this);
  if (is.null(name)) {
    msg <- paste(msg, " Header 'name' is missing.", sep="");
  } else if (nchar(name) == 0) {
    msg <- paste(msg, " Header 'name' is an empty string.", sep="");
  }

  msg <- trim(msg);
  if (nchar(msg) == 0)
    return(NULL);

  switch(action,
    error = throw(msg),
    warning = warning(msg)
  )

  msg; 
})



#########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the bioassay set"
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
setMethodS3("getName", "BaseFileBioassaysetinfo", function(this, ...) {
  getHeader(this, fields="name");
})



############################################################################
# HISTORY: 
# 2005-06-19
# o Added all Rdoc comments.
# o Created.
############################################################################
