###########################################################################/**
# @RdocClass PlatformInterface
#
# @title "The PlatformInterface class"
#
# \description{
#  @classhierarchy
#
#  A PlatformInterface provides methods for a given platform, e.g.
#  Affymetrix, Agilent, Illumina.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "Interface".}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("PlatformInterface", function(...) {
  extend(Interface(), "PlatformInterface");
})


setMethodS3("byName", "PlatformInterface", function(static, name, ...) {
  className <- sprintf("%sPlatform", capitalize(name));
  clazz <- Class$forName(className);
  newInstance(clazz);
}, static=TRUE);


setMethodS3("getName", "PlatformInterface", function(this, ...) {
  name <- grep("Platform", class(this), value=TRUE);
  name <- name[1];
  name <- gsub("Platform", "", name);
  name;
})

setMethodS3("findUnitNamesFile", "PlatformInterface", abstract=TRUE);

setMethodS3("getUnitNamesFile", "PlatformInterface", abstract=TRUE);



setMethodS3("getAromaUgpFile", "PlatformInterface", function(static, ...) {
  AromaUgpFile$byName(...);
}, static=TRUE)




############################################################################
# HISTORY:
# 2008-05-18
# o Created.
############################################################################
