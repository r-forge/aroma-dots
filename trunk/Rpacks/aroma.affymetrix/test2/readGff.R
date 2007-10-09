setMethodS3("readGff", "default", function(file, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  parseData <- function(data, ...) {
    # Data fields: 
    #  <seqname> <source> <feature> <start> <end> 
    #  <score> <strand> <frame> [attributes] [comments]
    colClasses <- c(seqname="character", source="character", feature="character", start="integer", end="integer", score="double", strand="character", frame="integer", attributes="character", comments="character"); 

    con <- textConnection(data, local=TRUE);
    on.exit(close(con));

    data <- read.table(con, colClasses="character", header=FALSE, sep="\t", fill=TRUE);
    colnames(data) <- names(colClasses)[1:ncol(data)];

    # Coerce integers
    for (mode in c("integer", "double")) {
      for (cc in which(colClasses == mode)) {
        values <- data[,cc];
        values <- trim(values);
        values[values == "."] <- NA;
        storage.mode(values) <- mode;
        data[,cc] <- values;
      }
    }

    data;    
  } # parseData()



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'file':
  if (is.character(file)) {
    con <- file(file, open="r");
    on.exit(close(con));
  } else if (inherits(file, "connection")) {
    con <- file;
  } else {
    stop("Argument 'file' must be a pathname or a connection: ", class(file)[1]);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Parse file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  tracks <- list();
  track <- NULL;
  data <- NULL;
  while (TRUE) {
    line <- readLines(con=con, n=1);
    if (length(line) == 0)
      break;

    # Strip comments
    line <- gsub("^#.*$", "", line);


    # A special header?
    if (regexpr("^browser", line) != -1) {
      # Ignore
    } else if (regexpr("^track", line) != -1) {
      # Store old track?
      if (!is.null(track)) {
        track$data <- parseData(data);
        key <- track$attributes$name;
        tracks[[key]] <- track;
      }

      tail <- gsub("^track", "", line);
      tail <- trim(tail);
      pattern <- "([a-zA-Z0-9_]*)=(\"[^\"]*\"|[^ ]*)(.*)";
      attrs <- list();
      while (nchar(tail) > 0) {
        key <- gsub(pattern, "\\1", tail);
        value <- gsub(pattern, "\\2", tail);
        value <- gsub("^\"", "", value);
        value <- gsub("\"$", "", value);
        attrs[[key]] <- value;
        tail <- gsub(pattern, "\\3", tail);
        tail <- trim(tail);
      }

      # A new track
      track <- list(attributes=attrs);
      data <- NULL;
    } else {
      data <- c(data, line);
    }
  }

  # Store last track
  if (!is.null(track)) {
    track$data <- parseData(data);
    key <- track$attributes$name;
    tracks[[key]] <- track;
  }

  tracks;  
}) # readGff()


############################################################################
# HISTORY:
# 2007-10-08
# o Created.
############################################################################
