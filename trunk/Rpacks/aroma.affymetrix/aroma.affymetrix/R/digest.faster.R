digest <- function(object, algo=c("md5", "sha1", "crc32"), serialize=TRUE, file=FALSE, length=Inf) {
  algo <- match.arg(algo);
  if (is.infinite(length))
    length <- -1;

  if (serialize && !file) {
    object <- serialize(object, connection=NULL, ascii=TRUE);
#    object <- paste(object, collapse="");
    object <- rawToChar(object);
  } else if (!is.character(object)) {
    stop("Argument object must be of type character if serialize is FALSE");
  }
  object <- as.character(object);
  algoint <- switch(algo, md5=1, sha1=2, crc32=3);
  if (file) {
    algoint <- algoint + 100;
    if (file.access(object) < 0) {
      stop(c("Can't open input file:", object));
    }
  }

  .Call("digest", as.character(object), as.integer(algoint),
                                      as.integer(length), PACKAGE="digest");
} # digest()


############################################################################
# HISTORY:
# 2007-01-06
# o Created.
############################################################################
