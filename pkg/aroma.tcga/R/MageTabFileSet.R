setConstructorS3("MageTabFileSet", function(...) {
  extend(TabularTextFileSet(...), "MageTabFileSet");
})


setMethodS3("byPath", "MageTabFileSet", function(static, ..., pattern=NULL, fileClass=getFileClass(static)) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fileClass':
  clazz <- Class$forName(fileClass);
  dfStatic <- getStaticInstance(clazz);
  dfStatic <- Arguments$getInstanceOf(dfStatic, getFileClass(static), .name="fileClass");

  # Argument 'pattern':
  if (!is.null(pattern)) {
    pattern <- Arguments$getRegularExpression(pattern);
  }

  # Default filename extension pattern
  if (is.null(pattern)) {
    pattern <- getExtensionPattern(dfStatic);
  }

  byPath.GenericDataFileSet(static, ..., pattern=pattern, fileClass=fileClass);
}, static=TRUE);



############################################################################
# HISTORY:
# 2009-10-25
# o Created.
############################################################################
