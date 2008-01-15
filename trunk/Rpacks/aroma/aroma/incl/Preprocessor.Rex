setConstructorS3("Preprocessor", function(mad=NULL, method=NULL, args=NULL) {
  if (!is.null(mad)) {
    if (!inherits(mad, "MicroarrayData"))
      throw("Object to be preprocessed is not of class MicroarrayData: ",
            getClass(mad));
    if (!is.null(method)) {
      if (!is.character(method))
        throw("Argument 'method' must be  : ",
            getClass(mad));
  }
  extend(Object(), "Preprocessor",
    .mad = mad
  )
})


setMethodS3("as.character", "Preprocessor", function(this) {
  s <- paste(getClass(this), 
             ": MicroarrayData object to be preprocessed is: ",
             as.character(this$.mad), sep="");
  s;
})








############################################################################
# HISTORY:
# 2002-10-27
# o Created!
############################################################################
