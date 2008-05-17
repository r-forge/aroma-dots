setMethodS3("writeColumnsToFiles", "TabularTextFile", function(this, destPath, filenameFmt="%s.txt", tags=NULL, columnName=NULL, header=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  writeHeaderComments0 <- function(con, hdr, commentPrefix="# ", ...) {
    hdr <- c(list(nbrOfHeaderRows=length(hdr)+1, hdr));    
    hdrStr <- unlist(hdr);
    hdrStr <- paste(names(hdrStr), hdrStr, sep="\t");
    hdrStr <- paste(commentPrefix, hdrStr, sep="");
    writeLines(con=con, hdrStr);
  }

  escapeFilename <- function(filename, ...) {
    filename <- gsub(":", "%3A", filename);
    filename <- gsub(";", "%3B", filename);
    filename;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'destPath':
  destPath <- Arguments$getWritablePath(destPath);

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- unlist(strsplit(tags, split=","));
    tags <- trim(tags);
    tags <- tags[nchar(tags) > 0];
  }

  # Argument 'filenameFmt':
  filenameFmt <- Arguments$getCharacter(filenameFmt);

  # Argument 'columnName':
  if (!is.null(columnName))
    columnName <- Arguments$getCharacter(columnName);

  # Argument 'header':
  if (is.null(header)) {
    header <- list(
      sourceFile=getFilename(this)
    );
  } else {
    header <- as.list(header);
  }

  hdrColumnName <- columnName;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify column names
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  columnNames <- getColumnNames(this);
  verbose && printf(verbose, "Column names [%d]:\n", length(columnNames));
  verbose && print(verbose, columnNames);
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract each column
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  colClassPatterns <- "character";
  for (cc in seq(along=columnNames)) {
    columnName <- columnNames[cc];
    verbose && enter(verbose, sprintf("Column #%d ('%s') of %d", 
                                       cc, columnName, length(columnNames)));
  
    fullname <- paste(c(columnName, tags), collapse=",");
    filename <- sprintf(filenameFmt, fullname);
    filename <- escapeFilename(filename);
    pathname <- file.path(destPath, filename);
    if (!isFile(pathname)) {
      names(colClassPatterns) <- sprintf("^%s$", columnName);
      values <- readDataFrame(this, colClassPatterns=colClassPatterns);
      values <- trim(values[[1]]);
      df <- data.frame(dummy=values, stringsAsFactors=FALSE);
      if (is.null(hdrColumnName)) {
        colnames(df) <- columnName;
      } else {
        colnames(df) <- hdrColumnName;
      }
      verbose && str(verbose, df);

      con <- file(pathname, open="w");
      header$createdOn <- format(Sys.time(), "%Y-%m-%d %H:%M:%S");
      header$column <- cc;
      header$columnName <- columnName;
      header$nbrOfDataRows <- nrow(df);
      writeHeaderComments0(con=con, header);
      write.table(file=con, df, quote=FALSE, sep="\t", row.names=FALSE);
      close(con);

      # Validate?
##      dbT <- GenericTabularFile(pathname);
##      print(dbT);
    } else {
      verbose && cat(verbose, "Column already extracted");
    }
  
    verbose && exit(verbose);
  } # for (cc ...)

  invisible(destPath);  
})


############################################################################
# HISTORY:
# 2008-05-05
# o Now some non-valid filename characters are escaped.
# o Added internal escapeFilename().
# 2008-05-01
# o Created.
############################################################################
