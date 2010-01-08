setMethodS3("readBCRSamplesTable", "TcgaDccDownloader", function(static, path="tcgaData", url="http://tcga-data.nci.nih.gov/tcga/BCRSamples.htm", ..., force=FALSE, verbose=FALSE) {
  require("XML") || throw("Package not loaded: XML");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getWritablePath(path);

  # Argument 'url':
  url <- Arguments$getCharacter(url);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving the BCR Biospecimen Barcode table");

  filename <- basename(url);
  fullname <- gsub("[.][^.]+$", "", filename);

  filenameHTML <- basename(url);
  pathnameHTML <- Arguments$getWritablePathname(filenameHTML, path=path);

  filenameTXT <- sprintf("%s.txt", fullname);
  pathnameTXT <- Arguments$getWritablePathname(filenameTXT, path=path);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load existing TXT file?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isFile(pathnameTXT)) {
    verbose && enter(verbose, "Reading existing TXT file");
    db <- TabularTextFile(pathnameTXT);
    verbose && print(verbose, db);
    data <- readDataFrame(db, colClassPattern=c("*"="character"));
    verbose && str(verbose, data);
    verbose && exit(verbose);
    verbose && exit(verbose);
    return(data);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Download HTML data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (force || !isFile(pathnameHTML)) {
    verbose && enter(verbose, "Downloading HTML file");
    download.file(url, destfile=pathnameHTML);
    pathnameHTML <- Arguments$getReadablePathname(pathnameHTML);
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Parse HTML file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Parsing HTML into XML");
  doc <- htmlParse(pathnameHTML, error=function(...){});
  verbose && str(verbose, doc);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract HTML table
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting table from HTML");
  path <- "/html/body//div[@class='eXtremeTable']/table[@id='ec_table']";
  table <- getNodeSet(doc, path);
  table <- table[[1]];
  verbose && str(verbose, table);

  # Sanity check
  stopifnot(length(table) > 0);

  # Extract the column names
  node <- getNodeSet(table, "//td[@class='tableHeader']");
  # Sanity check
  stopifnot(length(node) > 0);
  colnames <- sapply(node, FUN=xmlValue);
  ncol <- length(colnames);
  verbose && cat(verbose, "Number of columns: ", ncol, "\n", sep="");
  # Sanity check
  stopifnot(ncol < 999);
  verbose && cat(verbose, "Column names:");
  verbose && print(verbose, colnames);

  # Safe column names?
  colnames <- toCamelCase(colnames);


  # Extract table rows
  rows <- getNodeSet(table, "//tbody/tr");
  # Sanity check
  stopifnot(length(rows) > 0);
  nrow <- length(rows);
  verbose && cat(verbose, "Number of rows: ", nrow, "\n", sep="");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract HTML table as a data frame
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the body of all rows
  data <- lapply(rows, FUN=function(row) {
    cols <- getNodeSet(row, "td");
    sapply(cols, FUN=xmlValue, recursive=TRUE);
  });

  # Sanity checks
  ns <- sapply(data, FUN=length);
  uns <- unique(ns);
  stopifnot(length(uns) == 1);
  stopifnot(uns == ncol);

  # Turn into a matrix
  data <- unlist(data, use.names=FALSE);
  data <- matrix(data, ncol=ncol, byrow=TRUE);
  colnames(data) <- colnames;

  # Sanity check
  stopifnot(nrow(data) == nrow);
  stopifnot(ncol(data) == ncol);

  # Turn into a data frame
  data <- as.data.frame(data, stringsAsFactors=FALSE);
  verbose && str(verbose, data);

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save as TXT file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Writing TXT file");
  pathnameT <- sprintf("%s.tmp", pathnameTXT);
  pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=TRUE);
  verbose && cat(verbose, "Writing to temporary file: ", pathnameT);

  lastModified(pathnameHTML);
  header <- c();
  header <- c(header, sprintf("Import URL: %s", url));
  header <- c(header, sprintf("Imported on: %s", lastModified(pathnameHTML)));
  header <- c(header, sprintf("Imported file size: %s", file.info(pathnameHTML)$size));
  header <- c(header, sprintf("Imported using: %s v%s", getName(aroma.tcga), getVersion(aroma.tcga)));
  header <- c(header, sprintf("Number of data rows: %d", nrow(data)));
  comments <- paste("# ", header, sep="");
  comments <- paste(comments, collapse="\n");
  cat(file=pathnameT, comments, "\n", sep="");

  suppressWarnings({
    write.table(data, file=pathnameT, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, append=TRUE);
  });

  file.rename(pathnameT, pathnameTXT);
  if (isFile(pathnameT) || !isFile(pathnameTXT)) {
    throw("Failed to rename temporary file: ", pathnameT, " -> ", pathnameTXT);
  }
  verbose && exit(verbose);

  verbose && exit(verbose);

  data;
}, static=TRUE)


############################################################################
# HISTORY:
# 2010-01-07
# o Created.
############################################################################
