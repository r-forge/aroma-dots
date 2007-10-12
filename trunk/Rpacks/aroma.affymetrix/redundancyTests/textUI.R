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


textSelectFile <- function(path=".", pattern="[^~]$", ...) {
  pathHistory <- c();

  while(!isFile(path)) {
    pathHistory <- c(pathHistory, path);
    path <- Arguments$getReadablePathname(path);

    paths <- list.files(pattern=pattern, path=path, full.names=TRUE);
    if (length(paths) > 0) {
      options <- gsub(".*/", "", paths);
      options <- gsub(".(lnk|LNK)$", "", options);
      names(options) <- seq(along=options);

      if (length(pathHistory) > 1)
        options <- c(options, "-"="<back>");
      options <- c(options, "q"="<quit>");
      ans <- textMenu(options, title=path);

      if (options[ans] == "<quit>") {
        return(NULL);
      } else if (options[ans] == "<back>") {
        path <- pathHistory[length(pathHistory)-1];
        pathHistory <- pathHistory[seq(length=length(pathHistory)-2)];
      } else {
        path <- paths[ans];
      }
    } else {
      path <- paths[1];
    }
    path <- Arguments$getReadablePathname(path);
  }

  if (!isFile(path))
    return(NULL);

  path;
} # textSelectFile()


############################################################################
# HISTORY: 
# 2007-10-09
# o Added textSelectFile().
# 2007-01-11
# o One year aniversary of aroma.affymetrix!
# o Now mergeStrands and combineAlleles is automagically inferred when
#   defining a new CnChipEffectSet.
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
