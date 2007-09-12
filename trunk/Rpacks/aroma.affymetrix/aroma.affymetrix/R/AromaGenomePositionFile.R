setConstructorS3("AromaGenomePositionFile", function(...) {
  this <- extend(AffymetrixFile(...), "AromaGenomePositionFile");

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
}, abstract=TRUE)

setMethodS3("getFilenameExtension", "AromaGenomePositionFile", abstract=TRUE, protected=TRUE);


setMethodS3("readHeader", "AromaGenomePositionFile", function(this, ...) {
  # Open file
  pathname <- getPathname(this);
  con <- file(pathname, open="rb");
  on.exit(close(con));

  # Number of bytes in file header
  H <- readBin(con=con, what=integer(), size=4, signed=FALSE);

  # Skip file header
  seek(con=con, where=H, origin="current", rw="r");

  # Number of elements
  J <- readBin(con=con, what=integer(), size=4, signed=FALSE);

  # Number of bytes per chromosome value
  C <- readBin(con=con, what=integer(), size=1, signed=FALSE);

    # Number of bytes per position value
  P <- readBin(con=con, what=integer(), size=1, signed=FALSE);

  # Validation
  if (J == 0)
    throw("Genome position file format error: Number of genome positions is zero.");

  if (!C %in% c(1,2)) {
    throw("Genome position file format error: Number of bytes per chromosome value is not in {1,2}: ", C);
  }

  if (!P %in% c(1,2,4)) {
    throw("Genome position file format error: Number of bytes per position value is not in {1,2,4}: ", P);
  }

 list(nbrOfBytesInHeader=H, nbrOfElements=J, bytesPerChromosome=C, bytesPerPosition=P);
}, private=TRUE)


setMethodS3("nbrOfElements", "AromaGenomePositionFile", function(this, ...) {
  hdr <- readHeader(this);
  hdr$nbrOfElements;
}, protected=TRUE)

setMethodS3("getChipType", "AromaGenomePositionFile", function(this, ...) {
  getName(this, ...);
})

setMethodS3("getCdf", "AromaGenomePositionFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    chipType <- getChipType(this);
    cdf <- AffymetrixCdfFile$fromChipType(chipType);
    this$.cdf <- cdf;
  }
  
  cdf;
})

setMethodS3("getGenomeVersion", "AromaGenomePositionFile", function(this, ...) {
  tags <- getTags(this, ...);
  tags <- grep("^hg", tags, value=TRUE);
  tags;
}, protected=TRUE)

setMethodS3("create", "AromaGenomePositionFile", function(static, chipType, tags=NULL, nbrOfElements, path=NULL, overwrite=FALSE, ...) {
  # Argument 'nbrOfElements':
  nbrOfElements <- Arguments$getInteger(nbrOfElements, range=c(0,100e6));

  # Argument 'path':
  path <- Arguments$getWritablePath(path);

  # Generate pathname
  fullname <- paste(c(chipType, tags), collapse=",");
  filename <- sprintf("%s.%s", fullname, getFilenameExtension(static));
  pathname <- Arguments$getWritablePathname(filename, path=path);

  if (overwrite || !isFile(pathname)) {
    # Number of bytes in file header (for now always an empty file header)
    H <- as.integer(0);
  
    J <- nbrOfElements;
    C <- as.integer(1);
    P <- as.integer(4);
    
    # Open file
    con <- file(pathname, open="wb");
    on.exit(close(con));
  
    # Write file header
    writeBin(con=con, H, size=4, endian="little");
    seek(con=con, where=H, origin="current", rw="w");
  
    # Write data header
    writeBin(con=con, J, size=4, endian="little");
    writeBin(con=con, C, size=1, endian="little");
    writeBin(con=con, P, size=1, endian="little");
  
    # Write empty data
    # Chromosome index: Write NAs
    naValue <- as.integer(256^C-1);
    value <- rep(naValue, J);
    writeBin(con=con, value, size=C, endian="little");
    
    # Position: Write NAs
    suppressWarnings({
      # For P == 4, this will return NA, but that's ok.
      naValue <- as.integer(256^P-1);
    })
    value <- rep(naValue, nbrOfElements);
    writeBin(con=con, value, size=P, endian="little");
  }

  # Return
  newInstance(static, pathname);
}, static=TRUE, protected=TRUE)



setMethodS3("readData", "AromaGenomePositionFile", function(this, idxs=NULL, fields=c("chromosome", "position"), ..., verbose=FALSE) {
  # Data header
  hdr <- readHeader(this);
  H <- hdr$nbrOfBytesInHeader;
  J <- hdr$nbrOfElements;
  C <- hdr$bytesPerChromosome;
  P <- hdr$bytesPerPosition;
  
  # Open file
  pathname <- getPathname(this);
  con <- file(pathname, open="rb");
  on.exit(close(con));

  dataOffset <- 10+H;

  # Argument 'idxs':
  if (is.null(idxs)) {
    idxs <- seq_len(J);
  } else if (is.logical(idxs)) {
    idxs <- which(idxs);
    idxs <- Arguments$getIndices(idxs, range=c(1, J));
  } else {
    idxs <- Arguments$getIndices(idxs, range=c(1, J));
  }

  colClasses <- rep("integer", length(fields));
  names(colClasses) <- fields;
  data <- dataFrame(colClasses=colClasses, nrow=length(idxs));
  rownames(data) <- make.unique(as.character(idxs));
  
  # Chromosome indices
  if ("chromosome" %in% fields) {
    offset <- dataOffset;
    seek(con=con, where=offset, rw="r");
    values <- readBin(con=con, n=J, what=integer(), size=C, signed=FALSE, endian="little");
    values <- values[idxs];
    naValue <- 256^C - 1;
    values[values == naValue] <- NA;
    data[,"chromosome"] <- values;
    rm(values);
  }
  
  # Positions
  if ("position" %in% fields) {
    offset <- dataOffset + J*C;
    seek(con=con, where=offset, rw="r");
    values <- readBin(con=con, n=J, what=integer(), size=P, signed=FALSE, endian="little");
    values <- values[idxs];
    naValue <- 256^P - 1;
    values[values == naValue] <- NA;
    data[,"position"] <- values;
    rm(values);
  }
  
  data;
}, protected=TRUE)


setMethodS3("updateData", "AromaGenomePositionFile", function(this, idxs=NULL, chromosome=NULL, position=NULL, ..., verbose=FALSE) {
  # Data header
  hdr <- readHeader(this);
  H <- hdr$nbrOfBytesInHeader;
  J <- hdr$nbrOfElements;
  C <- hdr$bytesPerChromosome;
  P <- hdr$bytesPerPosition;
  
  # Open file
  pathname <- getPathname(this);
  con <- file(pathname, open="r+b");
  on.exit(close(con));

  dataOffset <- 10+H;

  # Argument 'idxs':
  if (is.null(idxs)) {
    idxs <- seq_len(J);
  } else {
    idxs <- Arguments$getIndices(idxs, range=c(1, J));
  }

  # Argument' chromosome':
  if (!is.null(chromosome)) {
    chromosome <- Arguments$getIndices(chromosome, range=c(0, 256^C-1), disallow=c("NaN"));
    chromosome <- rep(chromosome, length.out=length(idxs));
  }
  
  # Argument' position':
  if (!is.null(position)) {
    position <- Arguments$getIndices(position, range=c(0, 256^P-1), disallow=c("NaN"));
    position <- rep(position, length.out=length(idxs));
  }
  
  
  # Chromosome indices
  if (!is.null(chromosome)) {
    offset <- dataOffset;
    seek(con=con, where=offset, origin="start", rw="r");
    values <- readBin(con=con, n=J, what=integer(), size=C, signed=FALSE, endian="little");
    values[idxs] <- chromosome;
    # Update NAs
    naValue <- 256^C - 1;
    values[is.na(values)] <- naValue;
    suppressWarnings({
      values <- as.integer(values);
    })
    # Write data
    seek(con=con, where=offset, origin="start", rw="w");
    writeBin(con=con, values, size=C, endian="little");
    rm(values);
  }
  
  # Positions
  if (!is.null(position)) {
    offset <- dataOffset + J*C;
    seek(con=con, where=offset, origin="start", rw="r");
    values <- readBin(con=con, n=J, what=integer(), size=P, signed=FALSE, endian="little");
    values[idxs] <- position;
    # Update NAs
    naValue <- 256^P - 1;
    values[is.na(values)] <- naValue;
    suppressWarnings({
      values <- as.integer(values);
    })
    # Write data
    seek(con=con, where=offset, origin="start", rw="w");
    writeBin(con=con, values, size=P, endian="little");
    rm(values);
  }

  invisible(this);
}, protected=TRUE)


setMethodS3("indexOfElements", "AromaGenomePositionFile", function(this, names, ...) {
  throw("Do not know how to look up element names for class ", class(this), ".");
}, protected=TRUE)


setMethodS3("[", "AromaGenomePositionFile", function(this, i=NULL, j=1:2, drop=FALSE) {
  # Argument 'i':
  if (!is.null(i)) {
    if (is.character(i)) {
      i <- indexOfElements(this, i);
    }
  }

  # Argument 'j':
  if (is.null(j)) {
    j <- 1:2;
  } else if (is.character(j)) {
    fields <- j;
  }
  
  if (is.numeric(j)) {
    fields <- c("chromosome", "position")[j];
  }
  
  # Read data
  data <- readData(this, idxs=i, fields=fields);

  # Drop dimensions?
  if (drop && length(dim(data)))
    data <- drop(data);
  
  data;
})


setMethodS3("subset", "AromaGenomePositionFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  data <- readData(this);
  subset(data, ...);
})

setMethodS3("getElementsAt", "AromaGenomePositionFile", function(this, chromosome, range=NULL, ..., verbose=FALSE) {
  # Stratify by chromosome
  data <- readData(this, fields="chromosome")[[1]];
  keep <- !is.na(data) & (data %in% chromosome);
  idxs <- which(keep);

  if (!is.null(range)) {
    data <- readData(this, idxs=idxs, fields="position")[[1]];
    keep <- !is.na(data);
    keep <- keep & (range[1] <= data & data <= range[2]);
    idxs <- idxs[keep];
  }
  
  idxs;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2007-09-11
# o Moved getCdf() from AromaUgpFile to AromaGenomePositionFile.
# 2007-03-04
# o Now readData() returns a data frame. 
# o Added subset().
# o Added file header, that is, reserved a slot in the file format for 
#   future support of a file header.
# 2007-03-03
# o Created. Can import genome information data.
############################################################################
