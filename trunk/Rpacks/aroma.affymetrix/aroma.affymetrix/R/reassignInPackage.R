###########################################################################/**
# @RdocDefault reassignInPackage
#
# @title "Re-assigns a new value to an existing object in a loaded package"
#
# \description{
#  @get "title".
# }
# 
# @synopsis 
#
# \arguments{
#   \item{name}{The name of the object to be replaced."}
#   \item{pkgName}{The name of the package where the object lives."}
#   \item{value}{The new value to be assigned.}
#   \item{keepOld}{If @TRUE, the old value is kept as an attribute in
#     the new object.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) the new object.
# }
#
# @author
#
# \seealso{
#   See \code{assignInNamespace()} in @see "utils::getFromNamespace".
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("reassignInPackage", "default", function(name, pkgName, value, keepOld=TRUE, ...) {
  # Get the environment where to look for the function to replace
  envName <- sprintf("package:%s", pkgName);
  if (!envName %in% search())
    throw("Package not loaded: ", pkgName);
  env <- as.environment(envName);

  # Argument 'name':    
  if (!exists(name, envir=env)) {
    throw("Argument 'name' does not refere to an existing object: ", name);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Patch
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the object to be replaced
  oldValue <- get(name, envir=env);

  # Get environment of this object
  env <- environment(oldValue);

  # Assign this environment to the new object
  environment(value) <- env;

  # Keep the old value?
  if (keepOld)
    attr(value, "oldValue") <- oldValue;

  unlockBinding(name, env);
  assignInNamespace(name, value, ns=pkgName, envir=env);
  assign(name, value, envir=env);
  lockBinding(name, env);

  invisible(value);
}, private=TRUE) # reassignInPackage()

############################################################################
# HISTORY:
# 2007-07-04
# o This function might end up in the R.utils package later.
# o Added Rdoc comments.
# o reassignInPackage() is a generic function to replace an old object in a
#   locked package namespace.
# o Created from .patchAffymetrix.R.
############################################################################  
