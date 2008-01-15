#########################################################################/**
# @RdocClass GalLayout
#
# @title "Class representing GenePix Array List (GAL) data"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments accepted by the Layout constructor.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# \examples{
#  \dontrun{
#   layout <- GalLayout$read("fish.gal")
#   genes <- getGeneGroups(layout)
#   print(genes)
#  }
#
#  \dontrun{
#   # Gives:
#    [1] "GeneGroups: 7769 groups defined. 7680 genes with 1 replicates, 75 
#    genes with 3 replicates, 5 genes with 6 replicates, 4 genes with 16 
#    replicates, 4 genes with 64 replicates, 1 genes (3XSSC) with 193 
#    replicates."
#  }
# }
#
# \seealso{
#  The superclass @see "Layout".
# }
#
# \references{
#  The GenePix Array List (GAL) file format,
#  \url{http://www.axon.com/GN_GenePix_File_Formats.html}.
# }
#*/#########################################################################
setConstructorS3("GalLayout", function(...) {
  extend(Layout(...), "GalLayout"
  )
})


#########################################################################/**
# @RdocMethod read
#
# @title "Reads a GenePix Array List (GAL) file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The name of the file.}
#   \item{path}{Optional path where the data should be written.}
#   \item{solve}{If @TRUE, the reading is slower but more forgiving to
#         errors in the input file.}
#   \item{verbose}{If @TRUE, helpful information is printed at each step
#         of the parsing of the input file.}
# }
#
# @author
#
# \examples{\dontrun{See help(GalLayout) for an example.}}
#
# \seealso{
#   For writing a GalLayout object back to a file @seemethod "write".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("read", "GalLayout", function(this, filename, path=NULL, solve=FALSE, verbose=FALSE) {
  filename <- Arguments$getReadablePathname(filename, path);  

  if (verbose) cat("Reading file ", filename, "...", sep="");
  
  fh <- file(filename, "r");

  #=========================================================================
  # GAL Header
  #=========================================================================

  #-------------------------------------------------------------------------
  # Line 1: "ATF   1.0"
  # File type and version number.
  #-------------------------------------------------------------------------
  lines <- scan(fh, what=list("",""), sep="\t", quote="", nlines=1, strip.white=TRUE, flush=TRUE, quiet=TRUE);
  fileFmt <- lines[[1]];
  if (fileFmt != "ATF") {
    warning(paste("Unknown file format: ", fileFmt, "\n", sep=""));
  }

  # Currently the ATF version number can be "1.0" or "1".
  version <- lines[[2]];
  if (version != "1.0" && version != "1") {
    warning(paste("Unknown ATF file format version: ", fileFmt, " ", version, "\n", sep=""));
  }

  if (verbose)
    cat(fileFmt, version, "\n");

  #-------------------------------------------------------------------------
  # Line 2: "24   43"
  # Number of optional header records and number of data fields (columns).
  #-------------------------------------------------------------------------
  lines <- scan(fh, what=list(0,0), sep="\t", quote="", nlines=1, strip.white=TRUE, flush=TRUE, quiet=TRUE);
  nbr.of.optional.headers <- lines[[1]];
  nbr.of.columns <- lines[[2]];

  if (verbose)
    cat(nbr.of.optional.headers, nbr.of.columns,"\n");

  if (verbose)
    cat("Reading", nbr.of.optional.headers, "optional header records:\n");
  settings <- c();
  for (k in seq(nbr.of.optional.headers)) {
    #-----------------------------------------------------------------------
    # Line 3:2+nbr.of.optional.headers: 
    # "FieldName=Field Value"
    #-----------------------------------------------------------------------
    lines <- scan(fh, what=list(""), sep="\t", quote="\"", nlines=1, strip.white=TRUE, flush=TRUE, quiet=TRUE);
    record <- strsplit(lines[[1]], "=");
    field.name <- record[[1]][1];
    field.value <- record[[1]][2];
    settings[field.name] <- field.value;
    if (verbose)
      cat(field.name, "=", field.value, "\n", sep="");
  }

  #-------------------------------------------------------------------------
  # Line 2+nbr.of.optional.headers+1:end
  # The actual data itself
  #-------------------------------------------------------------------------
  knownHeaders <- c("Block"="integer", "Column"="integer", "Row"="integer", "Name"="character", "ID"="character");

  # Prescan the headers to identify version
  if (verbose)
    cat("Reading field names...");
  # scan() can't read column names that contains non ASCII 0-127 characters!
  header <- readLines(con=fh, n=1);
  close(fh);
  
  # Remove leading blanks
#  header <- gsub("^[ ]*", "", header);  # Does not work. Why? /HB 020807
  # Remove trailing blanks
  header <- gsub("[ ]*$", "", header);
  # Remove quotation marks
  header <- gsub("\"", "", header);

  # Extract the column names
  header <- unlist(strsplit(header, "\t"));
  if (verbose)
    cat("ok\n");

  galVersion <- NA;
  unknownHeaders <- c();
  colClasses <- c();
  for (k in seq(along=header)) {
    colNames <- names(knownHeaders);
    typeIdx <- which(sapply(seq(colNames), FUN=function(i) regexpr(colNames[i], header[k]) == 1));
    if (length(typeIdx) == 0) {
      unknownHeaders <- c(unknownHeaders, k);
      colClasses <- c(colClasses, NA);
    } else {
      typeIdx <- typeIdx[1];  # Bullet proof!
      colClasses <- c(colClasses, knownHeaders[typeIdx]);
    }
  }

  if (length(unknownHeaders) > 0) {
    unknown <- paste(header[unknownHeaders], collapse=", ");
    warning(paste(sep="", "GenePix file format warning: Some of the field names are not recognized. It might be because it is a new/old version. The unknown fields are: ", unknown, "."));
  }

  if (verbose)
    cat("Reading result table...");
  # NOTE: read.table(fh,...) does NOT work! Probably a bug in the skip part!
  # /Henrik 2001-03-07
  if (solve == TRUE) {
    # It is not possible to set colClasses <- NULL, which gives an error.
    df <- read.table(filename, sep="\t", quote="", skip=nbr.of.optional.headers+3, header=FALSE, na=c("NA", "NaN", "Error"));
  } else {
    df <- read.table(filename, sep="\t", quote="", skip=nbr.of.optional.headers+3, colClasses=colClasses, header=FALSE, na=c("NA", "NaN", "Error"));
  }
  if (verbose)
    cat(nrow(df), "lines. ok!\n");
  colnames(df) <- header;

  # "Name" and "ID" becomes factors, which are "harder" to work with.
  # By changing them to vectors the size increases (just) a little bit.
  # Change all factor levels that should be character's.
  for (k in which(colClasses == "character"))
    df[[k]] <- as.character(df[[k]]);

  # Change all factor levels that should be double's.
  for (k in which(colClasses == "double"))
    df[[k]] <- as.numeric(as.character(df[[k]]));
  
  if (verbose)
    cat("Successfully read self file.\n");

  #-------------------------------------------------------------------------
  # Extracting additional field to make life easier...
  #-------------------------------------------------------------------------
  nbrOfGenes <- length(df[,1]);
  nbrOfGrids <- max(df$Block);
  nspot.r    <- max(df$Row);
  nspot.c    <- max(df$Column);
  
  #-------------------------------------------------------------------------
  # Figuring out the layout structure of the microarray
  #-------------------------------------------------------------------------
  # Tries to figure out the number of rows of grids from the x coordinates.

  # Number of spots in each grid
  gridSize <- nspot.r * nspot.c;

  # The (x,y) position of each grid can be found in the header variables
  # 'Block1', 'Block2', ..., 'BlockN' where N is the number of grids.
  # The column of each of the BlockX variables might be separated 
  # with ',' (commas) or '\t' (tab).
  blockVarNames <- paste("Block", 1:nbrOfGrids, sep="");
  xyPos <- NULL;
  for (name in blockVarNames) {
    value <- settings[name];
    values <- strsplit(value, "[,\t]")[[1]];
    if (nchar(values[1]) == 0)
      values <- values[-1];
    values <- as.numeric(values);
    xyPos <- rbind(xyPos, c(x=values[1], y=values[2]));
  }

  # Figure out where the grids changes row
  dx <- sign(sign(c(0,xyPos[,"x"])-c(xyPos[,"x"],0))+1);

  ngrid.r <- sum(dx);
  ngrid.c <- max(df$Block) / ngrid.r;

  # Assert that out calculation are correct.
  if (round(ngrid.c) != ngrid.c)
    warning("Number of grid rows seems to be wrong! Please report this error to henrikb@braju.com.\n");

  # Make sure the name's and the id's are not factors!
  name <- as.character(df[,"Name"]);
  id <- as.character(df[,"ID"]);

  this <- GalLayout(ngrid.r, ngrid.c, nspot.r, nspot.c, name=name, id=id);
  this$version  <- version;
  this$settings <- settings;

  this;
}, static=TRUE)  # read()


setMethodS3("getFieldNames", "GalLayout", function(this) {
  c("Block", "Row", "Column", "ID", "Name");
})

setMethodS3("nbrOfFields", "GalLayout", function(this) {
  length(getFieldNames(this));
})

setMethodS3("getBlock", "GalLayout", function(this) {
  rep(1:nbrOfGrids(this), each=gridSize(this))
})

setMethodS3("getGridRow", "GalLayout", function(this) {
  ((getBlock(this)-1) %/% this$ngrid.c) + 1;
})

setMethodS3("getGridColumn", "GalLayout", function(this) {
  ((getBlock(this)-1) %% this$ngrid.c) + 1;
})

setMethodS3("getRow", "GalLayout", function(this) {
  row <- matrix(1:this$nspot.r, nrow=this$nspot.r, ncol=this$nspot.c, byrow=FALSE);
  row <- as.vector(t(row));
  row <- rep(row, times=this$ngrid.r*this$ngrid.c);
  col <- matrix(1:this$nspot.c, nrow=this$nspot.c, ncol=this$nspot.r, byrow=FALSE);
  col <- as.vector(col);
  col <- rep(col, times=this$ngrid.r*this$ngrid.c);
  row;
})

setMethodS3("getColumn", "GalLayout", function(this) {
  row <- matrix(1:this$nspot.r, nrow=this$nspot.r, ncol=this$nspot.c, byrow=FALSE);
  row <- as.vector(t(row));
  row <- rep(row, times=this$ngrid.r*this$ngrid.c);
  col <- matrix(1:this$nspot.c, nrow=this$nspot.c, ncol=this$nspot.r, byrow=FALSE);
  col <- as.vector(col);
  col <- rep(col, times=this$ngrid.r*this$ngrid.c);
  col;
})

setMethodS3("as.data.frame", "GalLayout", function(x) {
  # To please R CMD check...
  this <- x;

  df <- NULL;
  fieldNames <- getFieldNames(this);
  for (k in seq(along=fieldNames)) {
    name <- fieldNames[k];
    if (k == 1)
      df <- data.frame(this[[name]])
    else
      df <- cbind(df, data.frame(this[[name]]));
  }
  colnames(df) <- fieldNames;
  df;
})


setMethodS3("write", "GalLayout", function(this, filename, path=NULL, overwrite=FALSE, quote=FALSE, verbose=FALSE) {
  filename <- Arguments$getWritablePathname(filename, path, mustNotExist=!overwrite);  

  knownHeaders <- c("Block"="integer", "Column"="integer", "Row"="integer", "Name"="character", "ID"="character");

  settings <- this$settings;
  if (!is.element("Creator", names(settings)))
    settings <- c(settings, Creator=paste("com.braju.sma, http://www.braju.com/R/,", date()));

  fh <- file(filename, "w");
  on.exit(close(fh)); # Will be called even on an error.

  # Write ATF header
  cat("ATF\t", this$version, "\n", file=fh);
  cat(file=fh, length(settings), "\t", nbrOfFields(this), "\n", sep="");

  # Write GAL headers
  settingsStr <- paste("\"", names(settings), "=", settings, "\"\n", sep="");
  cat(file=fh, settingsStr, append=TRUE, sep="");
  rm(settings, settingsStr);

  # Get the column classes for each data field.
  colNames <- getFieldNames(this);
  knownNames <- names(knownHeaders);
  colClasses <- rep(NA, length(colNames));
  for (k in 1:length(colNames)) {
    found <- FALSE;
    for (l in seq(knownNames)) {
      found <- (regexpr(knownNames[l], colNames[k]) != -1);
      if (found) break;
    }
    if (found)
      colClasses[k] <- knownHeaders[l];
  }
  rm(knownNames, knownHeaders);
  
  stringCols <- (!is.na(colClasses) & colClasses == "character");

  # Write GAL data header
  colNames <- paste(paste("\"", colNames, "\"", sep=""), collapse="\t");
  cat(colNames, "\n", file=fh, append=TRUE);
  rm(colNames);

  # Get the actual data.
  df <- as.data.frame(this);
  write.table(df, file=fh, append=TRUE, sep="\t", eol="\n", row.names=FALSE, col.names=FALSE, quote=quote);
})  # write()


############################################################################
# HISTORY:
# 2005-10-21
# o Replace 'overwrite' arguments with 'mustNotExist' in calls to Arguments. 
# 2005-07-19
# o Replaced all path="" arguments to path=NULL.
# 2005-06-11
# o Making use of Arguments in R.utils.
# 2003-09-19
# o Added Rdoc comments for read().
# 2002-12-11
# o BUG FIX: Update read() to read both comma AND tab delimited headers.
# 2002-11-20
# o Added read() and write(), which were tested and verified with the
#   fish.gal GAL file.
# o Created.
############################################################################

