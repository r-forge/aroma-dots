setConstructorS3("AromaPlatformInterface", function(...) {
  extend(Interface(), "AromaPlatformInterface");
})


setMethodS3("getPlatform", "AromaPlatformInterface", abstract=TRUE);


setMethodS3("getAromaPlatform", "AromaPlatformInterface", function(this, ..., force=FALSE) {
  ap <- this$.ap;

  if (force || is.null(ap)) {
    platform <- getPlatform(this, ...);
    ap <- AromaPlatform$byName(platform, ...);
    this$.ap <- ap;
  }

  ap;
})



############################################################################
# HISTORY:
# 2008-06-12
# o Added abstract getPlatform().
# 2008-05-18
# o Created.
############################################################################
