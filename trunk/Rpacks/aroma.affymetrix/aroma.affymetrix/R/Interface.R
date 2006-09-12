###########################################################################/**
# @RdocClass Interface
#
# @title "The Interface class"
#
# \description{
#  @classhierarchy
#
#  This class represents a special set of classes whose purpose is to
#  provide methods (but not fields) shared by multiple different classes.
# }
# 
# @synopsis
#
# \arguments{
#   \item{core}{The core value.}
#   \item{...}{Not used.}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("Interface", function(core=NA, ...) {
  this <- core;
  class(this) <- "Interface";
  this;
})


setMethodS3("extend", "Interface", function(this, ...className, ...) {
  class(this) <- unique(c(...className, class(this)));
  this;
})


setMethodS3("uses", "Interface", function(this, ...) {
  setdiff(class(this), "Interface");
})


setMethodS3("as.character", "Interface", function(this, ...) {
  # Check if there are class "after" this one
  pos <- which("Interface" == class(this));
  isLast <- (pos == length(class(this)));
  if (isLast) {
    s <- paste(class(this), collapse=", ");
  } else {
    s <- NextMethod("as.character", this, ...);
  }
  s;
})


setMethodS3("print", "Interface", function(this, ...) {
  print(as.character(this), ...);
})


############################################################################
# HISTORY:
# 2006-09-11
# o Added trial version of an Interface class.
############################################################################
