digest <- function(object, algo=c("md5", "sha1", "crc32"), serialize=TRUE, file=FALSE, length=Inf) {
  algo <- match.arg(algo);
  if (is.infinite(length))
    length <- -1;

  if (serialize && !file) {
    object <- serialize(object, connection=NULL, ascii=TRUE);

    # Exclude serialization header (non-data dependent bytes but R version
    # specific).  In ASCII, the header consists of for rows ending with
    # a newline ('\n').  We remove these.
    # The end of 4th row should be within the first 30 bytes (typically 18).
    nHeader <- which(object[1:30] == as.raw(10))[4];
    # Remove header
    object <- object[-(1:nHeader)];

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

.assertDigest <- function(onDiff=c("error", "warning", "message"), ...) {
  # Argument 'onDiff':
  onDiff <- match.arg(onDiff);

  # Test against the digest() patch v0.4.2
  d0 <- "78a10a7e5929f8c605f71823203c0dc5";
  d1 <- digest(0);
  if (!identical(d1, d0)) {
    msg <- sprintf("Assertion failed: Detected inconsistency in digest(0) (%s != %s). The effect of this is that the generated cache files will be named differently on this platform/R version than in another.", d1, d0);
    if (onDiff == "error") {
      throw(msg);
    } else if (onDiff == "warning") {
      warning(msg);
    } else {
      cat(msg);
    }
  }
} # .assertDigest()


############################################################################
# HISTORY:
# 2007-03-08
# o Thanks to Luke Tierney's reply on my r-devel question of the serialize
#   header, we now look for the 4th newline, which is more robust to do
#   when serializing to ASCII.
# 2007-03-07
# o Added .assertDigest().
# o BUG FIX: serialize() gives different results depending on R version.
#   The difference seems to be in raw bytes 8, 9 and 10.  In other words,
#   if those are excluded before passing the stream on to the "digester"
#   we get the same results.  From testing with "random" object we also
#   know that there are at most 18 bytes in the header.
# 2007-01-06
# o Created.
############################################################################
