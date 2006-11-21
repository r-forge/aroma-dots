# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Generic user interface
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
selectMenu <- function(choices, title="Select/unselect items (0 when done)", ...) {
  selected <- rep(FALSE, length(choices));
  ans <- -1;
  while (ans != 0) {
    currChoices <- paste(c("[ ]", "[x]")[as.integer(selected)+1], choices, sep=" ");
    ans <- menu(choices=currChoices, title=title, ...);
    if (ans > 0)
      selected[ans] <- !selected[ans];
  }
  choices[selected];
} # selectMenu()


selectOrder <- function(choices, title="Select items one by one (0 to keep rest)", ...) {
  res <- c();
  while (length(choices) > 1) {
    if (length(res) > 0) {
      msg <- paste(seq(along=res), ": ", res, sep="");
      msg <- paste(msg, collapse=", ");
      msg <- paste("Currently selected items: ", msg, "\n", title);
    } else {
      msg <- title;
    }
    ans <- menu(choices=choices, title=msg, ...);
    if (ans == 0)
      break;
    res <- c(res, choices[ans]);
    choices <- choices[-ans];
  }
  res <- c(res, choices);
  res;
} # selectOrder()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data set(s) for with chip type(s) to be combined in one.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
selectDataSets <- function(paths="raw", pattern=NULL, ...) {
  paths <- sapply(paths, FUN=filePath, expandLinks="any");
  paths <- list.files(pattern="_(199|200)[0-9]", path=paths, full.names=TRUE);
  paths <- paths[sapply(paths, FUN=isDirectory)];
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Scan for chip types
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  paths <- sapply(paths, FUN=function(path) {
    pattern <- "^(Mapping|BI_SNP)"
    list.files(pattern=pattern, path=path);
  })
  keep <- sapply(paths, FUN=function(x) length(x) > 0);
  paths <- paths[keep];
  print(paths);
  
  uChipTypes <- sort(unique(unlist(paths)));
  print(uChipTypes);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Select chip type
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  chipTypes <- selectMenu(uChipTypes);
  
  # Filter out data sets
  keep <- sapply(paths, FUN=function(x) all(chipTypes %in% x));
  paths <- names(paths)[keep];
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Select data sets
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  paths <- selectMenu(paths);
  paths <- selectOrder(paths);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Define data sets
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  names <- NULL;
  name <- NA;
  tags <- NULL;
  dataSets <- list();
  for (kk in seq(along=chipTypes)) {
    chipType <- chipTypes[[kk]];
    files <- list();
    for (path in paths) {
      path <- filePath(path, chipType, expandLinks="any");
      if (!isDirectory(path))
        next;
  
      print(path);
      ds <- AffymetrixCelSet$fromFiles(path);
      names <- c(names, getName(ds));

      files <- append(files, as.list(ds));
    }
    ds <- AffymetrixCelSet(files=files);

    # When joining several data sets, we have to get a new name.
    if (length(names) > 0) {
      while (identical(name, NA)) {
        ans <- menu(choices=c("<new name>", names), title="Several data sets were joined. Choose the name you want to use for the new data set.");
        if (ans == 1) {
          while(TRUE) {
            name <- readline("Enter new name (blank for menu): ");
            name <- trim(name);
            if (name == "") {
              name <- NA;
              break;
            } else {
              tryCatch({
                setName(ds, name);
                break;
              }, error=function(ex) {});
            }
          } # while(TRUE)
        } else {
          name <- names[ans-1];
        }
      }
      setName(ds, name);

      # Give the option to add tags to joined data sets
      if (kk == 1) {
        tags <- readline("Enter additional tags (separated by commas or space): ");
        tags <- trim(tags);
        tags <- unlist(strsplit(tags, split="[ ,]"));
        tags <- trim(tags);
        tags <- unique(tags);
        if (tags == "")
          tags <- NULL;
      }
    }

    if (!is.null(tags))
      setTags(ds, c(getTags(ds), tags));

    dataSets[[chipType]] <- ds;
  } # for (chipType ...)

  dataSets;
} # selectDataSets()

