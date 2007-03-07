setMethodS3("getAttributes", "AffymetrixFile", function(this, ...) {
  attrs <- this$.attributes;
  if (is.null(attrs)) {
    attrs <- list();
  } else {
    # Always return attributes in lexicographic order by names
    o <- order(names(attrs));
    attrs <- attrs[o];
  }
  attrs;
}, protected=TRUE)

setMethodS3("setAttributes", "AffymetrixFile", function(this, ...) {
  # Argument '...':
  args <- list(...);
  names <- names(args);
  if (is.null(names)) {
    throw("No named arguments specified.");
  }
  
  # Update the attributes.
  attrs <- this$.attributes;
  attrs[names] <- args;
  this$.attributes <- attrs;

  invisible(args);
}, protected=TRUE)

setMethodS3("getAttribute", "AffymetrixFile", function(this, name, defaultValue=NULL, ...) {
  attrs <- this$.attributes;
  if (name %in% names(attrs)) {
    value <- attrs[[name]];
  } else {
    value <- defaultValue;
  }
  value;
}, protected=TRUE)

setMethodS3("setAttribute", "AffymetrixFile", function(this, name, value, ...) {
  attrs <- this$.attributes;
  attrs[[name]] <- value;
  this$.attributes <- attrs;

  invisible(attrs[name]);
}, protected=TRUE)

setMethodS3("testAttributes", "AffymetrixFile", function(this, select, ...) {
  # Get the attributes to be tested
  attrs <- getAttributes(this);
  expr <- substitute(select);
  res <- eval(expr, envir=attrs, enclos=parent.frame());
  res;
}, protected=TRUE)

setMethodS3("setAttributesByTags", "AffymetrixFile", function(this, tags=getTags(this), ...) {
  newAttrs <- list();

  # Get all <name>=<value> tags
  pattern <- "^([A-z][A-z0-9]*)=(.*)$";
  values <- grep(pattern, tags, value=TRUE);
  for (kk in seq(along=values)) {
    tag <- values[[kk]];
    key <- gsub(pattern, "\\1", tag);
    value <- gsub(pattern, "\\2", tag);

    # Try to coerce:
    suppressWarnings({
      value2 <- as.integer(value);
      if (!identical(value2 == value, TRUE)) {
        value2 <- as.double(value);
        if (!identical(value2 == value, TRUE)) {
          value2 <- as.character(value);
        }
      }
      value <- value2;
    })

    newAttrs <- c(newAttrs, setAttribute(this, key, value));
  }

  # Return updated attributes
  invisible(newAttrs);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2007-03-05
# o Added parseTagsAsAttributes(), which now also tries to coerce values.
# o Added support for (in-memory) attributes.
# 2007-02-07
# o Added getChecksum(), writeChecksum(), readChecksum(), and 
#   compareChecksum() and validateChecksum(). I did this because I noticed 
#   by chance that some of my CEL files transferred via an external HDD got
#   corrupt probe signals.
# 2007-01-14
# o Added a test for "unknown" (=unused) arguments to constructor.
# 2007-01-07
# o Added hasTags() and hasTag().
# 2006-11-02
# o Added getFullName(), getTags() and redefined getName().
# 2006-09-15
# o Added stextSize().
# 2006-08-27
# o Added stextLabel() and stextLabels(). stext is for "side text", cf. 
#   mtext for "margin text". stext() is slightly more convenient than mtext
#   when it comes to different font sizes.
# o Added copyTo().
# 2006-08-14
# o Added abstract fromFile().
# 2006-08-11
# o Created from AffymetrixDataFile in order to represent CDF files too.
############################################################################
