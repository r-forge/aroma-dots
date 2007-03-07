digest <- function(object, algo=c("md5", "sha1", "crc32"), serialize=TRUE, file=FALSE, length=Inf) {
  algo <- match.arg(algo);
  if (is.infinite(length))
    length <- -1;

  if (serialize && !file) {
    object <- serialize(object, connection=NULL, ascii=TRUE);
    # Exclude serialize "header", because it looks like it contains
    # R version specific information.
    object <- object[-(1:10)];
    object <- rawToChar(object);
  } else if (!is.character(object)) {
    stop("Argument object must be of type character if serialize is FALSE");
  }
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
# 2007-03-07
# o BUG FIX: serialize() gives different results depending on R version.
#   The difference seems to be in raw bytes 8, 9 and 10.  In other words,
#   if those are excluded before passing the stream on to the "digester"
#   we get the same results.  To be on the safe side, we exclude the
#   first ten raw bytes.
# 2007-01-06
# o Created.
############################################################################
