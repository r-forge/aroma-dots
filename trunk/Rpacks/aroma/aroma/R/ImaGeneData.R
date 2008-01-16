########################################################################/**
# @RdocClass ImaGeneData
#
# @title "The ImaGeneData class"
#
# @synopsis
#
# \description{
#  @classhierarchy
#
#  Creates an empty ImaGeneData object.
# }
#
# \arguments{
#   \item{layout}{A @see "Layout" object specifying the spot layout of the
#    slides in this data set.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \references{
#   BioDiscovery's ImaGene software, 
#   \url{http://www.biodiscovery.com/imagene.asp}
# }
#
# @examples "../incl/ImaGeneData.Rex"
#
# \note{
#   The laser beam with wavelength 635nm is the red laser and excites the
#   Cy5 dye, and the one with wavelength 532nm is the green laser which
#   excites the Cy3 dye.
# }
#*/#########################################################################
setConstructorS3("ImaGeneData", function(layout=NULL) {
  extend(MicroarrayData(layout=layout), "ImaGeneData",
    version=c(),
    settings=list()
  )
})




#########################################################################/**
# @RdocMethod read
#
# @title "Reads an ImaGene result file"
#
# \description{
#  Reads an ImaGene result file.
#  Currently ImaGene v4.1, v4.2, v5.0 files are supported. If the version
#  is not recognized it will try to make a best guess and read it in anyway.
#  However, it can not be guaranteed that all fields will be of the correct
#  data type.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file to be read.}
#   \item{path}{The path to the file.}
#   \item{verbose}{If @TRUE, information will printed out during
#                  the reading of the file.}
# }
#
# \value{Returns an @see "ImaGeneData" object.}
#
# @author
#
# \examples{
#   igG <- ImaGeneData$read("imagene234.cy3", path=system.file("data-ex", package="aroma"))
#   igR <- ImaGeneData$read("imagene234.cy5", path=system.file("data-ex", package="aroma"))
#   raw <- getRawData(igR, igG)
# }
#

# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("read", "ImaGeneData", function(this, filename, path=NULL, verbose=FALSE) {
  # Inner function parseHeader()
  parseHeader <- function(rows) {
    # Inner function parseSubtree()
    parseSubtree <- function(rows, currentRow=1) {
      isBeginTag <- function(row) {
        res <- (regexpr("^Begin[ ]+.*", row) != -1);
        if (res == TRUE) {
          row <- gsub("^Begin[ ]+", "", row);
          row <- gsub("[ \t]+$", "", row);
          attr(res, "type") <- row;
        }
        res;
      }
      
      isEndTag <- function(row) {
        res <- (regexpr("^End[ ]+.*", row) != -1);
        if (res == TRUE) {
          row <- gsub("^End[ ]+", "", row);
          row <- gsub("[ \t]+$", "", row);
          attr(res, "type") <- row;
        }
        res;
      }
     
      # Verify that first row is a Begin tag.
      row <- rows[currentRow];
      if (!isBeginTag(row))
        throw("Unexpected row contents: ", row);
      
      type <- attr(isBeginTag(row), "type")
    #  cat("Entering subtree ", type, "\n", sep="");
      
      tree <- list();
      while (currentRow < length(rows)) {
        currentRow <- currentRow + 1;
        row <- rows[currentRow];
    
        if (isBeginTag(row)) {
          # Entering a subtree
          subType <- attr(isBeginTag(row), "type")
          subTree <- parseSubtree(rows, currentRow=currentRow);
          currentRow <- attr(subTree, "endRow");
          tree[[subType]] <- subTree;
          rm(subTree);
        } else if (isEndTag(row)) {
          # Leaving a subtree
          endType <- attr(isEndTag(row), "type")
          if (endType != type)
            throw("ImaGene file format error: End tag (", endType, ") does not match begin tag (", type, ").");
    #      cat("Leaving subtree ", type, "\n", sep="");
          attr(tree, "endRow") <- currentRow;
          return(tree);
        } else {
          # Storing a field
          tree <- c(tree, row);
        }
      } # while(...)
    } # parseSubtree()
              
    tagToDataframe <- function(tag) {
      tag <- as.character(tag);   # 2004-01-23
      tag <- strsplit(tag, "\t");
      x <- matrix(unlist(tag[-1]), ncol=length(tag[[1]]), byrow=TRUE);
      colnames(x) <- tag[[1]];
      as.data.frame(x);
    }
    
    header <- parseSubtree(rows);
  
    # Fix all attributes
    for (k in seq(header)) {
      row <- header[[k]];
      if (!is.list(row)) {
        row <- as.character(row);
        row <- strsplit(row, "\t")[[1]];
        name <- row[1];
        value <- row[-1];
        header[[k]] <- value;
        names(header)[k] <- name;
      }
    }

    # Fix 'Field Dimensions' tag (not used for now)
    tag <- header[["Field Dimensions"]];
    if (!is.null(tag)) {
      df <- tagToDataframe(tag);
      for (col in c("Metarows", "Metacols", "Rows", "Cols"))
        df[[col]] <- as.numeric(as.character(df[[col]]));
      header[["Field Dimensions"]] <- df;
      rm(df, tag);
    }
    
    # Fix 'Measurement parameters' tag (not used for now)
    tag <- header[["Measurement parameters"]];
    if (!is.null(tag)) {
      df <- tagToDataframe(tag);
      for (col in seq(ncol(df)))
        df[[col]] <- as.character(df[[col]]);
      header[["Measurement parameters"]] <- df;
      rm(df, tag);
    }

    header;
  } # parseHeader()

  # Inner function readHeader()
  readHeader <- function(filename, beginRow, endRow, sep="\n", quote="") {
    fh <- file(filename, "r");
    on.exit(close(fh));
    skip <- beginRow-1;
    nlines <- endRow-skip;
    rows <- scan(fh, what=character(0), skip=skip, nlines=nlines, sep=sep, quote=quote, quiet=TRUE);
    # Remove all whitespace in the beginning of the rows
    rows <- gsub("^[ \t]+", "", rows);
    # Remove all whitespace at the end of the rows
    rows <- gsub("[ \t]+$", "", rows);
    header <- parseHeader(rows);
    header;
  } # readHeader()

  # Inner function readRawData()
  readRawData <- function(filename, beginRow, endRow, sep="\t", quote="") {
    field <- list("X"="character", "Field"="character", "Meta Row"="integer", "Meta Column"="integer", "Row"="integer", "Column"="integer", "Gene ID"="character", "Flag"="integer", "Signal Mean"="double", "Background Mean"="double", "Signal Median"="double", "Background Median"="double", "Signal Mode"="double", "Background Mode"="double", "Signal Area"="double", "Background Area"="double", "Signal Total"="double", "Background Total"="double", "Signal Stdev"="double", "Background Stdev"="double", "Shape Regularity"="double", "Empty"="double", "X Coord"="double", "Y Coord"="double", "Diameter"="double", "NA"=NA);
    # Fields added in v5.0
    field <- c(field, "Ignored Area"="double", "Spot Area"="double", "Ignored Median"="double", "Area To Perimeter"="double", "XCoord"="double", "YCoord"="double", "Position offset"="double", "Offset X"="double", "Offset Y"="double", "Expected X"="double", "Expected Y"="double", "CM-X"="double", "CM-Y"="double", "CM Offset"="double", "CM Offset-X"="double", "CM Offset-Y"="double", "Min Diam"="double", "Max Diam"="double", "Control"="character", "Failed Control"="integer", "Background contamination present"="integer", "Signal contamination present"="integer", "Ignored % failed"="integer", "Open perimeter failed"="integer", "Shape regularity failed"="integer", "Perim-to-area failed"="integer", "Offset failed"="integer", "Empty spot"="integer", "Negative spot"="integer");

    # Read the header line
    skip <- beginRow;
    df <- read.table(filename, header=TRUE, skip=skip, nrows=1, sep=sep, quote=quote, check.names=FALSE, comment.char="");
    header <- names(df);
    header[1] <- "X";
    unknown <- which(is.na(match(header, names(field))));
    if (length(unknown) > 0) {
      warning(paste("Unknown data field(s) found in ImaGene file. Will try to make an intelligent guess about the data type: ", paste(header[unknown], collapse=", "), sep=""));
      header[unknown] <- "NA";
    }

    colClasses <- unlist(field[header]);
    nrows <- endRow-skip-2;
    df <- read.table(filename, colClasses=colClasses, header=TRUE, skip=skip, nrows=nrows, sep=sep, quote=quote, check.names=FALSE, comment.char="");
    df <- df[,3:ncol(df)];
    colnames(df) <- gsub("[.]", " ", colnames(df));
    df;
  } # readRawData()

  filename <- Arguments$getReadablePathname(filename, path);  

  if (verbose) cat("Reading file ", filename, "...", sep="");

  # Support gzip'ed files too.
  if (regexpr("[.]gz$", filename) != -1) {
    tmpname <- tempfile();
    n <- gunzip(filename, tmpname);
    filename <- tmpname;
    on.exit(file.remove(tmpname));
  }
  
  fh <- file(filename, "r");

  ########################################################################
  # Step 1 - Prescan the file for Begin and End tags
  ########################################################################
  # Read the whole file into a number of strings containing the rows.
  rows <- readLines(fh);
  close(fh);

  # Find all Begin tags
  pos <- regexpr("^[ \t]*Begin[ ]*", rows);
  beginRow <- which(pos != -1);
  rows[beginRow] <- gsub("^[ \t]+", "", rows[beginRow]);
  rows[beginRow] <- gsub("[ \t]+$", "", rows[beginRow]);
  attr(beginRow, "match.length") <- attr(pos, "match.length")[beginRow];
  
  # Find all End tags
  pos <- regexpr("^[ \t]*End[ ]*", rows);
  endRow <- which(pos != -1);
  rows[endRow] <- gsub("^[ \t]+", "", rows[endRow]);
  rows[endRow] <- gsub("[ \t]+$", "", rows[endRow]);
  attr(endRow, "match.length") <- attr(pos, "match.length")[endRow];
  
  # Verify that the first Begin tag is a "Header"
  headerTag <- rows[beginRow[1]];
  if (headerTag != "Begin Header")
    throw("ImaGene file format error: Expected the first Begin tag to be a 'Begin Header' tag: ", headerTag);
  
  # Verify that the last End tag is a "End of File"
  len <- length(beginRow);
  endOfFileTag <- rows[endRow[len+1]];
  if (endOfFileTag != "End of File")
    throw("ImaGene file format error: The last End tag is not the End of File tag.");
  
  # Verify that the Begin tags match the End tags.
  if (len != length(endRow)-1)
    throw("ImaGene file format error: The number of Begin tags (", len, ") does not match the number of End tags (", length(endRow)-1, ").");
  
  beginTag <- substring(rows[beginRow], attr(beginRow, "match.length")+1);
  endTag <- substring(rows[endRow], attr(endRow, "match.length")+1)[-(len+1)];
  if (length(union(beginTag, endTag)) != len)
    throw("ImaGene file format error: The Begin tags do not match the End tags.");
  
  # Don't need the 'rows' anymore.
  rm(rows); # Save memory.

  ########################################################################
  # Step 2 - Read the 'Begin Header' tag
  ########################################################################
  tag <- "Header";
  begin <- beginRow[beginTag == tag];
  end   <- endRow[endTag == tag];
  header <- readHeader(filename, begin, end);

  # Verify that the Header tag a 'version' attribute.
  pos <- grep("^version$", names(header), ignore.case=TRUE)
  if (length(pos) == 0)
    throw("ImaGene file format error: Could not find 'version' attribute in the 'Header' tag.");
  version <- header[[pos]];
  if (!is.element(version, c("4.1", "4.2", "5.0")))
    warning(paste("Unsupported ImaGene file format version, but will give it a try anyway: ", version, sep=""));
  
  # Get the layout of the microarray
  fd <- header[["Field Dimensions"]];
  nbrOfSpots <- (fd$Metarows*fd$Metacols) * (fd$Rows*fd$Cols);
  
  ########################################################################
  # Step 3 - Read the 'Begin Raw Data' tag
  ########################################################################
  tag <- "Raw Data";
  begin <- beginRow[beginTag == tag];
  end   <- endRow[endTag == tag];
  nbrOfRows <- end-begin-2;
  if (nbrOfSpots != nbrOfRows) {
    warning("ImaGene file format error: According to the 'Field Dimensions' tag the cDNA microarray contains ", nbrOfSpots, " spots, but there are only ", nbrOfRows, " rows of data in the file. Will ignore the header information (anyway).");
  }
  
  if (verbose == TRUE)
    cat("Reading ", nbrOfRows, " rows of data...", sep="");
  df <- readRawData(filename, begin, end);
  if (verbose == TRUE)
    cat("done\n")


  ########################################################################
  # Step 4 - Creating the ImaGeneData object
  ########################################################################
  # Make sure the name's and the id's are not factors!
  id <- as.character(df[,"Gene ID"]);
  # Ignore the Header tag 'Field Dimensions'.
  ngrid.r <- max(df[,"Meta Row"]);
  ngrid.c <- max(df[,"Meta Column"]);
  nspot.r <- max(df[,"Row"]);
  nspot.c <- max(df[,"Column"]);
  layout <- Layout(ngrid.r,ngrid.c, nspot.r,nspot.c, id=id)
  rm(id);
  
  ig <- ImaGeneData(layout=layout)
  ig$version <- version;
  ig$header <- header;
  rm(header);

  ig$.fieldNames <- names(df);

  ig2 <- novirtual(ig)
  for (field in names(df)) {
    ig2[[field]] <- as.matrix(df[[field]]);
    df[[field]] <- NULL; # Save memory
  }

  rm(df); gc(); # To minimize memory usage!

  slidename <- basename(filename);
  setSlideNames(ig, slidename);

  ig;
}, static=TRUE) # read()


#########################################################################/**
# @RdocMethod write
#
# @title "Write an ImaGene result data file"
#
# \description{
#  Writes the ImaGeneData object to a file using ImaGene file format.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file to be written.}
#   \item{path}{The path to the file.}
#   \item{slide}{An @integer specifying which slide to be written to file.
#     Currently, only one slide at the time can be written.}
#   \item{...}{Arguments passed to \code{write.table}.}
# }
#
# \value{Returns nothing.}
#
# @author
#
# \examples{
#   ig <- ImaGeneData$read("imagene234.cy3", path=system.file("data-ex", package="aroma"))
#
#   # Writes the ImaGeneData object to a file named "temp.cy3" in a format
#   # that is as close as possible to the original format.
#   write(ig, "temp.cy3", overwrite=TRUE)
# }
#
# \seealso{
#   To read an ImaGene result data file see @seemethod "read".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("write", "ImaGeneData", function(this, filename, path=NULL, slide=1, overwrite=FALSE, ..., digits=13, verbose=FALSE) {
  filename <- Arguments$getWritablePathname(filename, path, mustNotExist=!overwrite);  

  slide <- validateArgumentSlide(this, slide=slide);

  writeTree <- function(tree, name, indent="", file=stdout()) {
    cat(file=file, indent, "Begin ", name, "\n", sep="", append=TRUE);
    indent <- paste(indent, "\t", sep="");
    names <- names(tree);
    for (k in seq(tree)) {
      nodeName <- names[k];
      node <- tree[[k]];
      if (is.list(node))
        writeTree(tree=node, name=nodeName, indent=indent, file=file)
      else
        cat(file=file, indent, node, "\n", sep="", append=TRUE);
    }
    indent <- substring(indent, 1, nchar(indent)-1);
    cat(file=file, indent, "End ", name, "\n", sep="", append=TRUE);
  } # writeTree()
  
  writeHeader <- function(header, file=stdout()) {
    cat(file=file, "Begin Header\n");
    indent <- "\t";
    headerNames <- names(header);
    for (k in seq(header)) {
      nodeName <- headerNames[k];
      node <- header[[k]];
      isList <- is.list(node);
      if (isList && nodeName == "Field Dimensions") {
        cat(file=file, indent, "Begin ", nodeName, "\n", sep="", append=TRUE);
        indent <- paste(indent, "\t", sep="");
        x <- paste(names(node), sep="\t", collapse="\t");
        cat(file=file, indent, x, "\n", sep="", append=TRUE);
        for (r in seq(nrow(node))) {
          cat(file=file, indent, append=TRUE);
          for (c in seq(ncol(node))) {
            cat(file=file, as.character(node[r,c]), "\t", append=TRUE);
          }
          cat(file=file, "\n", append=TRUE);
        }
        indent <- substring(indent, 1, nchar(indent)-1);
        cat(file=file, indent, "End ", nodeName, "\n", sep="", append=TRUE);
      } else if (isList && nodeName == "Measurement parameters") {
        cat(file=file, indent, "Begin ", nodeName, "\n", sep="", append=TRUE);
        indent <- paste(indent, "\t", sep="");
        x <- paste(names(node), sep="\t", collapse="\t");
        cat(file=file, indent, x, "\n", sep="", append=TRUE);
        for (r in seq(nrow(node))) {
          cat(file=file, indent, append=TRUE);
          for (c in seq(ncol(node))) {
            cat(file=file, as.character(node[r,c]), "\t", append=TRUE);
          }
          cat(file=file, "\n", append=TRUE);
        }
        indent <- substring(indent, 1, nchar(indent)-1);
        cat(file=file, indent, "End ", nodeName, "\n", sep="", append=TRUE);
      } else if (isList) {
        writeTree(tree=node, name=nodeName, indent=indent, file=file);
      } else {
        cat(file=file, indent, nodeName, "\t", node, "\n", sep="", append=TRUE);
      }
      rm(nodeName, node);
    }
    cat(file=file, "End Header\n", append=TRUE);
  } # writeHeader()

  writeRawData <- function(df, file=stdout()) {
    field <- list("X"="character", "Field"="character", "Meta Row"="integer", "Meta Column"="integer", "Row"="integer", "Column"="integer", "Gene ID"="character", "Flag"="integer", "Signal Mean"="double", "Background Mean"="double", "Signal Median"="double", "Background Median"="double", "Signal Mode"="double", "Background Mode"="double", "Signal Area"="double", "Background Area"="double", "Signal Total"="double", "Background Total"="double", "Signal Stdev"="double", "Background Stdev"="double", "Shape Regularity"="double", "Empty"="double", "X Coord"="double", "Y Coord"="double", "Diameter"="double", "NA"=NA, "indent"="character");
    # Fields added in v5.0
    field <- c(field, "Ignored Area"="double", "Spot Area"="double", "Ignored Median"="double", "Area To Perimeter"="double", "XCoord"="double", "YCoord"="double", "Position offset"="double", "Offset X"="double", "Offset Y"="double", "Expected X"="double", "Expected Y"="double", "CM-X"="double", "CM-Y"="double", "CM Offset"="double", "CM Offset-X"="double", "CM Offset-Y"="double", "Min Diam"="double", "Max Diam"="double", "Control"="integer", "Failed Control"="integer", "Background contamination present"="integer", "Signal contamination present"="integer", "Ignored % failed"="integer", "Open perimeter failed"="integer", "Shape regularity failed"="integer", "Perim-to-area failed"="integer", "Offset failed"="integer", "Empty spot"="integer", "Negative spot"="integer");
    
    header <- colnames(df);
    unknown <- which(is.na(match(header, names(field))));
    if (length(unknown) > 0) {
      warning(paste("Unknown data field(s) found in ImaGeneData object. Will try to make an intelligent guess about the data type: ", paste(header[unknown], collapse=", "), sep=""));
      header[unknown] <- "NA";
    }
    colClasses <- unlist(field[header]);
    
    cat(file=file, "Begin Raw Data\n", append=TRUE);
    colNames <- paste(header[-1], collapse="\t");
    cat(file=file, "\t", colNames, "\n", sep="", append=TRUE);
    write.table(df, file=file, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE, na="NaN", append=TRUE);
    cat(file=file, "End Raw Data\n", append=TRUE);
  } # writeRawData()

  fh <- file(filename, "w");
  on.exit(close(fh));

  # Write Header
  writeHeader(this$header, file=fh);

  # Write Raw Data
  df <- as.data.frame(this);
  indent <- rep("", nrow(df));
  field <- rep(as.character(this$header$Field$Field), nrow(df));
  df <- cbind(indent, Field=field, df[,3:ncol(df)]);
  rm(indent, field);
  oopt <- options(digits=digits);
  writeRawData(df, file=fh);
  options(oopt);
  rm(df);
  
  # Write End of File
  cat(file=fh, "End of File\n", append=TRUE);
}) # write()




#########################################################################/**
# @RdocMethod getRawData
#
# @title "Gets the raw intensites from two ImaGeneData objects"
#
# \description{
#  Extracts the red and green spot intensitites (both foreground and
#  background) from two ImaGeneData objects and returns a
#  @see "RawData" object.
# }
#
# @synopsis
#
# \arguments{
#   \item{igR, igG}{The "red" and the "green" ImageGeneData objects.}
#   \item{slides}{Specifying which slides to be extracted. If @NULL, 
#     all slides are considered.}
#   \item{fg}{If \code{"mean"}, the mean foreground intensities are returned.
#     If \code{"median"}, the median foreground intensities are returned.}
#   \item{bg}{If \code{"mean"}, the mean background intensities are returned.
#     If \code{"median"}, the median background intensities are returned.}
# }
#
# \value{
#   Returns a @see "RawData" object containing the specified slides.
# }
#
# \details{
#   The R and Rb channels will come from the igR object, and
#   the G and Gb channels will come from the igG object.
#   To swap the channels swap the arguments or use dyeSwap() afterwards.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getRawData", "ImaGeneData", function(igR, igG, slides=NULL, fg="median", bg="mean") {
  if (!inherits(igG, "ImaGeneData"))
    throw("Argument 'igG' is not of class ImaGeneData: ", data.class(igG));
  slides <- validateArgumentSlides(igR, slides=slides);
  slides <- validateArgumentSlides(igG, slides=slides);

  if (!equals(getLayout(igR), getLayout(igG)))
    throw("The two specified ImaGeneData objects (argument 'igR' and 'igG') do not have the same layout and can not be combined.");

  #-------------------------------------------------------------------------
  # RGs as in the sma library.
  #-------------------------------------------------------------------------
  if (fg == "mean")
    fgs <- c("Signal Mean")
  else if (fg == "median")
    fgs <- c("Signal Median")
  else
    throw("Argument 'fg' has an unknown value: ", fg);

  if (bg == "mean")
    bgs <- c("Background Mean")
  else if (bg == "median")
    bgs <- c("Background Median")
  else
    throw("Argument 'bg' has an unknown value: ", bg);

  nas <- which(is.na(match(c(fgs,bgs), getFieldNames(igR))));
  if (length(nas) > 0) {
    throw("Strange, some fields do not exists in the ImaGeneData object 'igR'. Please, report this error to the author of the package: ", c(fgs,bgs)[nas], collapse=", ");
  }
  nas <- which(is.na(match(c(fgs,bgs), getFieldNames(igG))));
  if (length(nas) > 0) {
    throw("Strange, some fields do not exists in the ImaGeneData object 'igR'. Please, report this error to the author of the package: ", c(fgs,bgs)[nas], collapse=", ");
  }

  R  <- igR[[fgs]][,slides];
  G  <- igG[[fgs]][,slides];
  Rb <- igR[[bgs]][,slides];
  Gb <- igG[[bgs]][,slides];
    
  RawData(R=R, G=G, Rb=Rb, Gb=Gb, layout=getLayout(igR), 
                                           extras=igR$.extras)
});

setMethodS3("as.RawData", "ImaGeneData", function(this, ...) {
  getRawData(this, ...);
})



setMethodS3("setLayout", "ImaGeneData", function(this, layout) {
  warning("For an ImaGeneData object it should not be necessary to set the layout explicitly. The layout is automatically calculated when the data is loaded from file.");
  NextMethod("setLayout", this);
})


setMethodS3("getArea", "ImaGeneData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include)) include <- seq(nbrOfSpots(this));
  as.matrix(this[["Signal Area"]][include,slides]);
})


setMethodS3("getBgArea", "ImaGeneData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include)) include <- seq(nbrOfSpots(this));
  as.matrix(this[["Background Area"]][include,slides]);
})


setMethodS3("getCircularity", "ImaGeneData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include)) include <- seq(nbrOfSpots(this));
  as.matrix(this[["Shape Regularity"]][include,slides]);
})


setMethodS3("getAbscent", "ImaGeneData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include)) include <- seq(nbrOfSpots(this));
  as.matrix((this[["Empty"]][include,slides] == 1));
})

setMethodS3("anonymize", "ImaGeneData", function(this, ...) {
  # Anonymize a copy of the Layout object; others might use this one.
  layout <- clone(getLayout(this));
  anonymize(layout, ...);
  # Assume the same layout on all slides!
  slides <- seq(this);
  this[["Gene ID"]][,slides] <- getID(layout);
})



setMethodS3("updateHeader", "ImaGeneData", function(this) {
  # Get the current header
#  layout <- c(Metarows=max(this[["Meta Row"]]), Metacols=max(this[["Meta Column"]]), Row=max(this[["Row"]]), Col=max(this[["Column"]]));
  layout <- getLayout(this);
  fd <- this$header[["Field Dimensions"]];
  fd[,"Metarows"] <- layout$ngrid.r;
  fd[,"Metacols"] <- layout$ngrid.c;
  fd[,"Rows"] <- layout$nspot.r;
  fd[,"Cols"] <- layout$nspot.c;
  this$header[["Field Dimensions"]] <- fd;
}, protected=TRUE)



#########################################################################/**
# @RdocMethod getSpotPosition
#
# @title "Gets physical positions of the spots"
#
# \description{
#  Gets physical positions (in pixels) of the spots on one or several
#  slides. 
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{Specifying which for slides the spot positions should
#     be extracted. If @NULL, all slides are considered.}
#   \item{index}{The spots for which the position is returned.
#     If @NULL all spots are considered.}
# }
#
# \value{Returns a @see "SpotPosition" object containing the
#   positions of the spots on the specified slides.}
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getSpotPosition", "ImaGeneData", function(this, slides=NULL, index=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (any(index < 1) || any(index > nbrOfSpots(this)))
    throw("Argument 'index' is out of range.");
  if (is.null(index))
    index <- 1:nbrOfSpots(this);

  xField <- "X Coord";
  yField <- "Y Coord";
  if (!hasField(this, xField) || !hasField(this, yField)) {
    throw("This ", data.class(this), " object is missing the fields '", xField, "' and/or '", yField, "', which specifies the physical position of the spots.");
  }

  x <- this[[xField]][index,slides];
  y <- this[[yField]][index,slides];
  SpotPosition(x=x, y=y);
})



setMethodS3("plotSpatial", "ImaGeneData", function(this, what=NULL, slide=1, pch=20, yaxt=NULL, col="auto", palette="redgreen", A.range=c(0,16), M.range=c(-1,1), style=c("real", "classic"), ...) {
  style = match.arg(style);
  if (style == "classic")
    return(plotSpatial.MicroarrayData(this, what=what, slide=slide, col=col, ...));
  
  if (length(slide) != 1 || !is.numeric(slide))
    throw("Argument 'slide' must be an integer.");
  
  if (!is.null(what) && !is.element(what, getFieldNames(this)))
    throw("Argument 'what' is refering to an unknown field: ", what);
  
  if (is.null(col) || col == "auto") {
    if (is.null(what))
      what <- "Signal Median";
    col <- getColors(this, what, slide=slide);
  }
  
  xy <- getSpotPosition(this, slides=slide);
  # Upper *left* corner is (0,0)
  x <- xy$x; y <- -xy$y; # Will also be used as xlab and ylab.
  
  if (missing(yaxt))
    yaxt0 <- "n";

  plot(x,y, pch=pch, yaxt=yaxt0, col=col, ...);
  
  if (missing(yaxt)) {
    # Add the y-axis with the correct ticks
    yaxp <- par("yaxp");
    yaxis <- seq(yaxp[1],yaxp[2], length=yaxp[3]+1);
    axis(2, at=yaxis, labels=-yaxis);
  }
})


setMethodS3("plotSpatial3d", "ImaGeneData", function(this, field="Signal Median", ...) {
  NextMethod("plotSpatial3d", this, field=field, ...);
})


############################################################################
# HISTORY:
# 2006-07-04
# o BUG FIX: read() of ImaGeneData could not handle '#' in the Gene ID.
#   Added comment.char="" to all read.table() calls.
# 2005-10-21
# o Replace 'overwrite' arguments with 'mustNotExist' in calls to Arguments. 
# 2005-07-19
# o Replaced all path="" arguments to path=NULL.
# 2005-06-11
# o Making use of Arguments in R.utils.
# 2004-08-15
# o Renamed as.RawData() to getRawData().
# 2004-03-09
# o BUG FIX: write() was mistakenly accepting more than once slide at the
#   time. It was a typo and now an error is thrown is more slides are given.
# 2004-01-23
# o BUG FIX: read() was crashing at a strsplit() call in the internal
#   tagToDataframe() function. Thanks to Brett Milash at the Huntsman 
#   Cancer Institute, US for providing the bug fix.
# 2003-12-04
# o Added a read-multiple-slides example to the ImaGeneData help. This was
#   done on request from Giovanni Martinelli, Italy.
# 2003-10-30
# o Added getSpotPosition().
# 2003-10-26
# o Added support to read *.gz files.
# 2003-10-19
# o BUG FIX: A typo made as.RawData() for ImaGeneData to return the *mean*
#   instead of the median.
# 2003-08-29
# o Added match.arg() where applicable.
# 2002-12-24
# o Added plotSpatial3d().
# 2002-12-21
# o When reading slides from files each slide is now named as the filename.
# 2002-12-05
# o BUG FIX: The get*() methods did not return a matrix if there was only
#   one slide.
# 2002-09-24
# o Changed the attribute 'path="."' to 'path=""' in read().
# o Added Rdoc comments for as.RawData().
# 2002-08-09
# o BUG FIX: A minor bug, but critical bug was found in the internal
#   function parseHeader() of the read() method. It assumed the hard coded
#   value "auto" in Header -> Measurement parameters -> Segmentation Method.
#   Modified the function to accept any kind of value since it is not used
#   anyway. Sangsoo, South Korea, reported the problem with "fixed circle".
# 2002-07-02
# o Can not write ImaGene v5.0 files also. Added an argument 'digits=13' to
#   write() to be more compatible with the output format ImaGene is
#   producing. Added an internal writeTree() function to write().
# o A fix in extract() in the MicroarrayData class asserts that the colnames
#   of a data frame are not renamed to 'safe names'. This fix made it
#   possible to keep the true column names when reading the data.
# o Finally I got the reply from Biodiscovery about how to download ImaGene.
#   In their demo version there are a few ImaGene v5.0 files, which I now
#   succesfully can load. Updated internal readRawData() in read().
# 2002-07-01
# o Spent about 10 hours to implement and test this class.
# o Added protected updateHeader().
# o Added support for write() and added its Rdoc comments. Verified to work
#   for both v4.1 and v4.2.
# o Added getArea(), getBgArea(), getCircularity(), getAbscent().
# o Added setLayout() and anonymize().
# o Added as.RawData() and its Rdoc comments.
# o Created the class ImaGeneData and its Rdoc comments. Can read
#   version 4.1 and 4.2. Noticed on the ImaGene webpage that v5.0 is out.
# 2002-06-07(?)
# o Register online to get a trial version of ImaGene and said that I wanted
#   it so I could get some information about the file format. Never got it!
#   I don't think Biodiscovery wants others to be able to read their files.
# 2002-05-29
# o Created. Did some initial trials to read ImaGeneData files.
# 2002-05-25
# o Asked support$biodiscovery.com for information about the ImaGene file
#   format, but I never got a reply.
############################################################################
