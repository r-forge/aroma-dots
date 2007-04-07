setMethodS3("readCfnUnits", "default", function(pathname, snps=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathname':
  pathname <- Arguments$getReadablePathname(pathname, mustExist=TRUE);

  # Argument 'units':
  if (!is.null(snps))
    snps <- Arguments$getIndices(snps);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Reading data from CFN file");
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve file header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hdr <- readCfnHeader(pathname);
  verbose && cat(verbose, "Header:");
  verbose && str(verbose, hdr);
  nbrOfBytes <- hdr$nbrOfBytes;
  nbrOfSnps <- hdr$nbrOfSnps;
  bytesPerSnp <- hdr$bytesPerSnp;

  # Validating units
  if (!is.null(snps))
    snps <- Arguments$getIndices(snps, range=c(1, nbrOfSnps));

  map <- matrix(1:(bytesPerSnp*nbrOfSnps), nrow=bytesPerSnp);
  map <- map + hdr$dataOffset;

  # Read subset of SNPs?
  if (!is.null(snps))
    map <- map[,snps,drop=FALSE];
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read all data
  raw <- readBin(pathname, what="raw", n=nbrOfBytes);

  # Bytes 13:16 contains (M) as floats
  rrM <- 13:16;
  # Bytes 17:20 contains (C) as integers
  rrC <- 17:20;
  rr <- setdiff(5:12, c(rrM, rrC));
#  rr <- c(4:19);
rr <- rrM;
  ncol <- length(rrM) / 4;
  map <- map[rrM,,drop=FALSE];
  theta <- readBin(raw[map], what="double", size=4, endian="little", 
                                                           n=ncol*ncol(map));
#  theta <- readBin(raw[map], what="integer", size=4, endian="little", 
  theta <- matrix(theta, ncol=ncol, byrow=TRUE);

#  colnames(theta) <- c("A", "B");

  verbose && exit(verbose);
  
  theta;
})


############################################################################
# HISTORY:
# 2007-04-06
# o Created.
############################################################################
