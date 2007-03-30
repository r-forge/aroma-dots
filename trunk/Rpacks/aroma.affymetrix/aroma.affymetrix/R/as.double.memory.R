setMethodS3("as.double", "double", function(x, ...) {
 if (is.null(attributes(x))) {
   x;
 } else {
   NextMethod("as.double", x, ...);
 }
})

setMethodS3("as.integer", "integer", function(x, ...) {
 if (is.null(attributes(x))) {
   x;
 } else {
   NextMethod("as.integer", x, ...);
 }
})

setMethodS3("as.raw", "raw", function(x, ...) {
 if (is.null(attributes(x))) {
   x;
 } else {
   NextMethod("as.raw", x, ...);
 }
})

setMethodS3("as.complex", "complex", function(x, ...) {
 if (is.null(attributes(x))) {
   x;
 } else {
   NextMethod("as.complex", x, ...);
 }
})

setMethodS3("as.logical", "logical", function(x, ...) {
 if (is.null(attributes(x))) {
   x;
 } else {
   NextMethod("as.logical", x, ...);
 }
})

setMethodS3("as.character", "character", function(x, ...) {
 if (is.null(attributes(x))) {
   x;
 } else {
   NextMethod("as.character", x, ...);
 }
})


############################################################################
# HISTORY:
# 2007-03-28
# o By default as.double() copies the object in memory, even if already
#   a double. This code avoids such waste of memory.
# o Created.
############################################################################
