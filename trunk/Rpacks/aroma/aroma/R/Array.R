# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Array
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
setConstructorS3("Array", function(array=NULL, dim=NULL) {
  this <- extend(BasicObject(array), "Array");
  if (is.null(dim))
    dim <- dim(array);

  # Tries to preserved the dimnames if the are compatible
  dimnames <- dimnames(array);
  dim(this) <- dim;
  compatible <- TRUE;
  for (k in seq(along=dim)) {
    names <- dimnames[[k]];
    if (!is.null(names) && length(names) != dim[k])
      compatible <- FALSE;
  }
  if (compatible)
    dimnames(this) <- dimnames;
  this;
})

setMethodS3("as.data.frame", "Array", function(x, ...) {
  # To please R CMD check...
  this <- x;

  as.data.frame(as.matrix(this), ...);
})


############################################################################
# HISTORY:
# 2002-11-27
# o Added as.data.frame().
# o Constructor now preserved the dimnames if they are compatible with the
#   dimensions.
# 2002-11-03
# o Created!
############################################################################
