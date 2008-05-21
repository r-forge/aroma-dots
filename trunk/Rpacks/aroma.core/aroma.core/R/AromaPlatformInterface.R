setConstructorS3("AromaPlatformInterface", function(...) {
  extend(Interface(), "AromaPlatformInterface");
})


setMethodS3("getAromaPlatform", "AromaPlatformInterface", function(this, ..., force=FALSE) {
  ap <- this$.ap;

  if (force || is.null(ap)) {
    platform <- getPlatform(this, ...);
    ap <- AromaPlatform$byName(platform, ...);
    this$.ap <- ap;
  }

  ap;
})



###########################################################################/**
# @RdocClass AromaPlatform
#
# @title "The AromaPlatform class"
#
# \description{
#  @classhierarchy
#
#  A AromaPlatform provides methods for a given platform, e.g.
#  Affymetrix, Agilent, Illumina.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("AromaPlatform", function(...) {
  extend(Object(), "AromaPlatform");
})


setMethodS3("byName", "AromaPlatform", function(static, name, ...) {
  className <- sprintf("%sPlatform", capitalize(name));
  clazz <- Class$forName(className);
  newInstance(clazz);
}, static=TRUE);


setMethodS3("getName", "AromaPlatform", function(this, ...) {
  name <- class(this)[1];
  name <- name[1];
  name <- gsub("Platform", "", name);
  name;
})

setMethodS3("findUnitNamesFile", "AromaPlatform", abstract=TRUE);

setMethodS3("getUnitNamesFile", "AromaPlatform", abstract=TRUE);


setMethodS3("getAromaUgpFile", "AromaPlatform", function(static, ...) {
  AromaUgpFile$byName(...);
}, static=TRUE)




############################################################################
# HISTORY:
# 2008-05-18
# o Created.
############################################################################
