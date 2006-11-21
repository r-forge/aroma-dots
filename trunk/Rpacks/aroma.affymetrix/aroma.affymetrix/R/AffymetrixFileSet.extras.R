setMethodS3("getDescription", "AffymetrixFileSet", function(this, ...) {
  path <- getPath(this);
  res <- list();
  for (kk in 1:3) {
    pathname <- file.path(path, "DESCRIPTION");
    if (isFile(pathname)) {
      tmp <- read.dcf(pathname);
      tmp <- as.list(as.data.frame(tmp));
      tmp <- lapply(tmp, FUN=as.character);
      for (kk in seq(along=tmp)) {
        key <- names(tmp)[kk];
        # Already assigned?
        if (key %in% names(res))
          next;
        res[[key]] <- tmp[[key]];
      }
      break;
    }
    path <- dirname(path);
  }

  res;
}, protected=TRUE)

setMethodS3("getIdentifier", "AffymetrixFileSet", function(this, ...) {
  path <- getPath(this);
  res <- NULL;
  for (kk in 1:3) {
    pathname <- file.path(path, "IDENTIFIER");
    if (isFile(pathname)) {
      res <- readLines(pathname);
      # Remove comments
      res <- trim(gsub("#.*", "", trim(res)));
      # Remove empty lines
      res <- res[nchar(res) > 0];
      break;
    }
    path <- dirname(path);
  }

  if (!is.null(res)) {
    res <- digest(list(res));
  }

  res;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2006-11-20
# o Made getDescription() protected.
# 2006-10-30
# o Added getDescription() which search and parse all DESCRIPTION files in
#   the data-set directory tree.
# o Added getIdentifier() which returns a 32-character long hexadecimal
#   hashcode for the "Identifier" string returned by getDescription().
#   If no such string exists, NULL is returned.  This will allow users
#   to specify their own identifiers.
############################################################################
