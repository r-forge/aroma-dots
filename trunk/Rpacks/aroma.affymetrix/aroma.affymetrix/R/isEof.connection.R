setMethodS3("isEof", "connection", function(con, ...) {
  # Remember position
  offset <- seek(con, rw="read");
  # Try to read next byte
  bfr <- readChar(con, nchars=1);
  # Reposition
  seek(con, where=offset, rw="read");
  # No more bytes?
  (nchar(bfr) == 0);
}) # isEof()


