# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Public functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
updateGraphics <- function(sets, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  hasPrefix <- function(name, prefix, ...) {
    (substring(name, 1, nchar(prefix)) == prefix);
  }

  hasAsterisk <- function(name, ...) {
    (regexpr("*", name, fixed=TRUE) != -1);
  }

  hasPlus <- function(name, ...) {
    (regexpr("+", name, fixed=TRUE) != -1);
  }

  hasHash <- function(name, ...) {
    (regexpr("#", name, fixed=TRUE) != -1);
  }

  for (kk in seq(along=sets)) {
    key <- names(sets)[kk];
    set <- sets[[kk]];
    name <- set$name;
    verbose && enter(verbose, sprintf("%d %s\n", kk, name));
  
    # Default colors
    col <- "black";
    lty <- 3;

    if (hasPrefix(name, "CRMA")) {
      col <- colors["CRMA"];
      if (hasPlus(name) & hasAsterisk(name) & hasHash(name)) {
        col <- "blue";
        lty <- 4;
      } else if (hasPlus(name) & !hasAsterisk(name)) {
        lty <- 1;
      } else if (hasPlus(name) & hasAsterisk(name)) {
        col <- "green";
        lty <- 2;
      } else if (hasHash(name) & !hasAsterisk(name)) {
        lty <- 3;
      } else {
#        col <- colors["APT"];
        lty <- 5;
      }
    } else if (hasPrefix(name, "dChip")) {
      col <- colors["dChip"];
      lty <- 1;
    } else if (hasPrefix(name, "APT")) {
      col <- colors["APT"];
      lty <- 1;
    } else if (hasPrefix(name, "GTC")) {
      col <- colors["GTC"];
      lty <- 1;
    }

    if (!identical(col, set$col)) {
      set$col <- col;
    }
    if (!identical(lty, set$lty)) {
      set$lty <- lty;
    }

    sets[[key]] <- set;

    rm(set, col, lty, name);
    verbose && exit(verbose);
  } # for (kk ...)

  sets;
} # updateGraphics()
