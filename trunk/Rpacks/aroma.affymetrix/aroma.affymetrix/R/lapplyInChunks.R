setMethodS3("lapplyInChunks", "numeric", function(idxs, fcn, ..., chunkSize=1, useNames=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fcn':
  if (!is.function(fcn)) {
    throw("Argument 'fcn' is not a function: ", mode(fcn));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "lapplyInChunks()");

  n <- length(idxs);

  # Nothing todo?
  if (length(n) == 0)
    return(list());
   
  remaining <- 1:n;
  chunk <- 1:chunkSize;
  nbrOfChunks <- ceiling(n / chunkSize);

  # Allocate return structure
  res <- vector("list", n);

  verbose && cat(verbose, "Number elements per chunk: ", chunkSize);

  count <- 1;
  while (length(remaining) > 0) {
    verbose && enter(verbose, "Chunk #", count, " of ", nbrOfChunks);
    # Last chunk?
    if (length(remaining) < chunkSize)
      chunk <- 1:length(remaining);

    # Indices for this chunk
    ii <- remaining[chunk];

    verbose && cat(verbose, "Elements: ");
    verbose && str(verbose, idxs[ii]);

    # Apply the function
    resChunk <- fcn(idxs[ii], ...);
    res[ii] <- resChunk;

    # Next chunk
    remaining <- remaining[-chunk];
    count <- count + 1;

    verbose && exit(verbose);
  } # while(length(idxs) > 0)

  # Add names?
  if (useNames)
    names(res) <- names(idxs);

  verbose && exit(verbose);

  res;
}) # lapplyInChunks()


############################################################################
# HISTORY:
# 2007-02-12
# o Created. 
############################################################################
