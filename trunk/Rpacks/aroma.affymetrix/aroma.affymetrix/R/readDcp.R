readDcp <- function(con, nbrOfUnits, fields=c("rawIntensities", "normalizedIntensities", "calls", "thetas", "thetaStds", "excludes"),  cells=NULL, units=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readElements <- function(con, idxs, nbrOfElements, size=1, skip=FALSE, ...) {
    # Local functions
    readBin2 <- function(con, what, size, signed, n) {
      if (what == "raw") {
        n <- n * size;
        size <- 1;
      }
#print(list(what=what, size=size, signed=signed, n=n));
      readBin(con, what=what, size=size, signed=signed, n=n);
    } # readBin2()


    if (skip) {
      seek(con, where=size*nbrOfElements, origin="current", rw="read")
      return(NULL);
    }

    # Read all elements?
    if (is.null(idxs)) {
      readBin2(con=con, size=size, ..., n=nbrOfElements);
    } else {
      r <- range(idxs);
      nbrOfElementsToSkip <- r[1]-1;
      nbrOfElementsToRead <- r[2]-r[1]+1;

      # Skip to first element to read?
      if (nbrOfElementsToSkip > 0) {
        seek(con=con, where=size*nbrOfElementsToSkip, 
                                           origin="current", rw="read");
      }

      # Read values
      values <- readBin2(con=con, size=size, ..., n=nbrOfElementsToRead);
      nbrOfBytesLeft <- size*(nbrOfElements - length(values));
      if (nbrOfBytesLeft > 0) {
        seek(con=con, where=nbrOfBytesLeft, origin="current", rw="read");
      }
      values;
    }
  } # readElements()


  readCells <- function(con, cells, what="integer", size=2, signed=FALSE, ...) {
    readElements(con=con, idxs=cells, nbrOfElements=nbrOfCells, 
                                   what=what, size=size, signed=signed, ...);
  }


  readUnits <- function(con, units, what="float", size=2, signed=FALSE, ...) {
    readElements(con=con, idxs=units, nbrOfElements=nbrOfUnits, 
                                   what=what, size=size, signed=signed, ...);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'con':
  if (is.character(con)) {
    pathname <- con;
    pathname <- Arguments$getReadablePathname(pathname, mustExist=TRUE);
    con <- file(con, open="rb");
    on.exit({
      if (!is.null(con))
        close(con);
      con <- NULL;
    });
  }

  if (!inherits(con, "connection")) {
    throw("Argument 'con' must be either a connection or a pathname: ", 
                                                            mode(con));
  }

  # Argument 'fields':  
  fields <- match.arg(fields, several.ok=TRUE);

  # Argument 'cells':
  if (!is.null(cells))
    cells <- Arguments$getIndices(cells);

  # Argument 'units':
  if (!is.null(units))
    units <- Arguments$getIndices(units);
  

  res <- list();

  res$header <- readDcpHeader(con=con, verbose=verbose);
#str(res);

  nbrOfCells <- res$header$CellDim^2;

  # Argument 'nbrOfUnits':
  nbrOfUnits <- Arguments$getInteger(nbrOfUnits, range=c(1, nbrOfCells));


  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, range=c(1, nbrOfCells));
  }

  if (!is.null(units)) {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits));
  }

  # Read intensities?
  for (field in c("rawIntensities", "normalizedIntensities")) {
    res[[field]] <- readCells(con, cells=cells, skip=(!field %in% fields));
#    str(res);
  }

  # Read calls?
  field <- "calls";
  res[[field]] <- readUnits(con, units=units, what="raw", size=1, 
                                               skip=(!field %in% fields));
#  str(res);

  # Read MBEI estimates?
  skip <- !any(c("thetas", "thetaStds", "excludes") %in% fields);
  raw <- readUnits(con, units=units, what="raw", size=12, skip=skip);
  if (!skip) {
    raw <- matrix(raw, nrow=12);
    n <- ncol(raw);
  
    field <- "thetas";
    if (field %in% fields) {
      res[[field]] <- readBin(con=raw[1:4,], what="double", size=4, signed=TRUE, n=n);
    }
    raw <- raw[-(1:4),,drop=FALSE];
  
    field <- "thetaStds";
    if (field %in% fields) {
      res[[field]] <- readBin(con=raw[1:4,], what="double", size=4, signed=TRUE, n=n);
    }
    raw <- raw[-(1:4),,drop=FALSE];
  
    field <- "excludes";
    if (field %in% fields) {
      res[[field]] <- readBin(con=raw[1:4,], what="integer", size=4, signed=TRUE, n=n);
    }
  }
  rm(raw);

# print(seek(con, read="r"));
#print(file.info(pathname2)$size);

  res;
} # readDcp()


##############################################################################
# HISTORY:
# 2008-01-30
# o Created.
##############################################################################

