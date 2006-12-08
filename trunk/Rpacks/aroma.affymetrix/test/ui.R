# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Generic user interface
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
textMenu <- function(choices, title="", prompt="Selection: ") {
  nc <- length(choices);
  keys <- names(choices);
  if (is.null(keys))
    keys <- seq_len(nc);
  op <- paste(format(keys, justify="right"), ": ", choices, sep="");

  if (nc > 10) {
    fop <- format(op);
    nw <- nchar(fop[1], "w") + 2;
    ncol <- getOption("width") %/% nw;
    if (ncol > 1) {
      op <- paste(fop, c(rep("  ", ncol-1), "\n"), sep="", collapse="");
    }
  }

  if (length(title) && nchar(title[1]))
    cat(title[1], "\n");
  cat("", op, "", sep="\n");

  keys <- trim(keys);
  repeat{
    ans <- readline(prompt);
    ans <- trim(ans);
    idx <- pmatch(ans, keys);
    if (is.finite(idx))
      return(idx);
    cat(gettext("Enter an item from the menu.\n"))
  }
} # textMenu()


selectMenu <- function(choices, selected=NULL, title="Select/unselect items", options=c("a!"="Select all", "n!"="Unselect all", "t!"="Toggle all", "q!"="Done"), header="%s (0 when done)", ...) {
  if (is.null(selected)) {
    selected <- rep(FALSE, length=length(choices));
  } else if (is.logical(selected)) {
    selected <- rep(selected, length.out=length(choices));
  } else if (is.numeric(selected)) {
    idx <- selected;
    selected <- rep(FALSE, length=length(choices));
    selected[idx] <- TRUE;
  } else if (is.character(selected)) {
    selected <- (choices %in% selected);
  }

  # Argument 'options':


  # Argument 'title' and 'header':
  title <- sprintf(header, title);

  nbrOfChoices <- length(choices);
  if (is.null(names(choices))) {
    names(choices) <- seq_len(nbrOfChoices);
  }


  repeat{
    currChoices <- paste(c("[ ]", "[x]")[as.integer(selected)+1], choices, sep=" ");
    names(currChoices) <- names(choices);
    optionIdxs <- nbrOfChoices + seq(along=options);
    ans <- textMenu(choices=c(currChoices, options), title=title, ...);
    if (ans > nbrOfChoices) {
      opt <- options[ans-nbrOfChoices];
      if (opt == "Select all") {
        selected[seq(along=currChoices)] <- TRUE;
      } else if (opt == "Unselect all") {
        selected[seq(along=currChoices)] <- FALSE;
      } else if (opt == "Toggle all") {
        selected <- !selected;
      } else if (opt == "Done") {
        break;
      }
    } else {
      selected[ans] <- !selected[ans];
    }
  }
  choices[selected];
} # selectMenu()



selectOrder <- function(choices, title="Select order", header="%s (0 to keep rest)", ...) {
  res <- c();
  while (length(choices) > 1) {
    if (length(res) > 0) {
      msg <- paste(seq(along=res), ": ", res, sep="");
      msg <- paste(msg, collapse=", ");
      msg <- paste("Currently selected items: ", msg, "\n");
      msg <- paste(msg, sprintf(header, title), "\n", sep="");
    } else {
      msg <- sprintf(header, title);
    }

    ans <- textMenu(choices=c(choices, "q!"="Done"), title=msg, ...);
    if (ans == length(choices)+1)
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
selectDataSets <- function(paths="raw", pattern=NULL, class=AffymetrixCelSet, ...) {
  paths <- sapply(paths, FUN=filePath, expandLinks="any");
  paths <- list.files(pattern=pattern, path=paths, full.names=TRUE);
  if (length(paths) == 0) {
    throw("Cannot select data set. No data sets found matching pattern '", 
                           pattern, "': ", paste(paths, collapse=", "));
  }
  paths <- paths[sapply(paths, FUN=isDirectory)];
  dataSetPaths <- paths;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Scan for chip types
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  paths <- lapply(paths, FUN=function(path) {
    pattern <- "^(Mapping|BI_SNP)"
    list.files(pattern=pattern, path=path);
  })
  names(paths) <- dataSetPaths;
  keep <- sapply(paths, FUN=function(x) length(x) > 0);
  paths <- paths[keep];
  
  uChipTypes <- sort(unique(unlist(paths)));
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Select chip type
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  chipTypes <- selectMenu(uChipTypes);
  if (length(chipTypes) == 0)
    chipTypes <- uChipTypes;

  # Filter out data sets
  keep <- sapply(paths, FUN=function(x) all(chipTypes %in% x));
  paths <- names(paths)[keep];
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Select data sets
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  names <- sort(basename(paths));
  names <- selectMenu(names);
  names <- selectOrder(names, title="Select order how data sets should be joined");
  idx <- match(names, basename(paths));
  paths <- paths[idx];
  
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
      ds <- class$fromFiles(path);
      names <- c(names, getName(ds));

      files <- append(files, as.list(ds));
    }
    ds <- class(files=files);

    # Remove duplicated arrays 
    if (!inherits(ds, "ChipEffectSet")) {
      dups <- isDuplicated(ds);
      if (any(dups)) {
        cat(sprintf("Removing %d duplicated arrays.\n", sum(dups)));
        ds <- extract(ds, !dups);
      }
    }

    # When joining several data sets, we have to get a new name.
    unames <- unique(names);
    if (length(unames) > 1) {
      while (identical(name, NA)) {
        ans <- textMenu(choices=c("<new name>", unames), title="Several data sets were joined. Choose the name you want to use for the new data set.");
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
          name <- unames[ans-1];
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
        if (identical(tags, ""))
          tags <- NULL;
      }
    }

    if (!is.null(tags))
      setTags(ds, c(getTags(ds), tags));

    dataSets[[chipType]] <- ds;
  } # for (chipType ...)

  dataSets;
} # selectDataSets()

############################################################################
# HISTORY: 
# 2006-12-02
# o Added textMenu().
# 2006-12-01
# o Now selectDataSets() uses only unique data sets names when asking for
#   a new name when merging several data sets, i.e. if there is only one
#   unique name, then that is used.
# o Now selectDataSets() removed duplicated arrays.
# 2006-11-27
# o Added argument 'selected' to selectMenu().
# 2006-11-22
# o Created.
############################################################################
