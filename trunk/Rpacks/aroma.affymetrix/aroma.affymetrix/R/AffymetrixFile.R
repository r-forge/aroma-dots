###########################################################################/**
# @RdocClass AffymetrixFile
#
# @title "The abstract AffymetrixFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixFile object represents a single Affymetrix file,
#  e.g. an Affymetrix CEL file or an Affymetrix CDF file.
#  Note that this class is abstract and can not be instanciated, but
#  instead you have to use one of the subclasses or the generic 
#  @seemethod "fromFile" method.
# }
# 
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file.}
#   \item{path}{An optional path to the file.}
#   \item{mustExist}{If @TRUE, an exception is thrown if the file does
#     not exists, otherwise not.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# \section{Filename convention}{
#   This package assumes a special filename structure allowing methods
#   to infer data-set names, sample names, file types, versions etc from 
#   the full path name.
#   
#   The filename of an \code{AffymetrixFile} is structured as follows:
#   \itemize{
#    \item{basename}{\code{"sample001,a,b,c.CEL"} 
#       (this follows the \R convention (but not the Unix convention)}
#    \item{name}{\code{"sample001"}}
#    \item{adjunct names}{\code{c("a", "b", "c")}}
#    \item{extension}{\code{"CEL"}}
#   }
# }
#
# @author
#
# \seealso{
#   An object of this class is typically part of an @see "AffymetrixFileSet".
# }
#*/###########################################################################
setConstructorS3("AffymetrixFile", function(filename=NULL, path=NULL, mustExist=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(filename)) {
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=mustExist);
  } else {
    pathname <- NULL;
  }

  extend(Object(), "AffymetrixFile",
    .pathname = pathname
  )
}, abstract=TRUE)


setMethodS3("clone", "AffymetrixFile", function(this, clear=TRUE, ...) {
  object <- NextMethod("clone", this, ...);
  if (clear)
    clearCache(object);
  object;
})


###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the Affymetrix file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("as.character", "AffymetrixFile", function(this, ...) {
  s <- paste(class(this)[1], ":", sep="");
  s <- paste(s, " Name: ", getName(this), ".", sep="");
  tags <- getTags(this);
  if (!is.null(tags)) {
    s <- paste(s, " Tags: ", paste(tags, collapse=","), ".", sep="");
  }
  s <- paste(s, " File type: ", getFileType(this), ".", sep="");
  s <- sprintf("%s Pathname: %s (%.2fMb).", s, getPathname(this), 
                                                   getFileSize(this)/1024^2);
  s <- sprintf("%s RAM: %.2fMb.", s, objectSize(this)/1024^2);
  s;
})


###########################################################################/**
# @RdocMethod getPathname
#
# @title "Gets the pathname of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPathname", "AffymetrixFile", function(this, ...) {
  this$.pathname;
})




###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path (directory) of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPath", "AffymetrixFile", function(this, ...) {
  dirname(this$.pathname);
})




###########################################################################/**
# @RdocMethod getFilename
#
# @title "Gets the filename of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character.
# }
#
# \details{
#  The filename of a file is the pathname excluding any path.
#  For instance, the filename of \code{path/to/foo,a.2,b.ext} is 
#  \code{foo,a.2,b.ext}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFilename", "AffymetrixFile", function(this, ...) {
  basename(this$.pathname);
})



###########################################################################/**
# @RdocMethod getFullName
#
# @title "Gets the full name of the file"
#
# \description{
#   @get "title", that is the filename without the extension.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character.
# }
#
# \details{
#  The full name of a file is the filename excluding any
#  extension (and period).
#  For instance, the full name of \code{path/to/foo,a.2,b.ext} is 
#  \code{foo,a.2,b}.
# }
#
# @author
#
# \seealso{
#   @seemethod "getName".
#   @seemethod "getTags".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFullName", "AffymetrixFile", function(this, ...) {
  name <- basename(this$.pathname);

  # Exclude filename extension
  name <- gsub("[.][a-zA-Z0-9][a-zA-Z0-9]*$", "", name);

  name;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character.
# }
#
# \details{
#  The name of a file is the part of the filename without the extension and
#  that preceeds any comma.
#  For instance, the name of \code{path/to/foo,a.2,b.ext} is \code{foo}.
# }
#
# @author
#
# \seealso{
#   @seemethod "getTags".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "AffymetrixFile", function(this, ...) {
  name <- getFullName(this, ...);

  # Keep anything before the first comma
  name <- gsub("[,].*$", "", name);
  
  name;
})




###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character @vector or @NULL.
# }
#
# \details{
#  The \emph{tags} of a filename are the comma separated parts of the
#  filename that follows the the first comma, if any, and that preceeds the
#  last period (the filename extension).
#  For instance, the tags of \code{path/to/foo,a.2,b.ext} are 
#  \code{a.2} and \code{b}.
# }
#
# @author
#
# \seealso{
#   @seemethod "getName".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "AffymetrixFile", function(this, ...) {
  name <- getFullName(this, ...);

  # Data-set name is anything before the first comma
  dsName <- gsub("[,].*$", "", name);

  # Keep anything after the data-set name (and the separator).
  name <- substring(name, nchar(dsName)+2);
  
  tags <- strsplit(name, split=",")[[1]];
  if (length(tags) == 0)
    tags <- NULL;
  tags;
})


setMethodS3("getLabel", "AffymetrixFile", function(this, ...) {
  label <- this$label;
  if (is.null(label))
    label <- getName(this, ...);
  label;
})

setMethodS3("setLabel", "AffymetrixFile", function(this, label, ...) {
  this$label <- label;
  invisible(this);
})



###########################################################################/**
# @RdocMethod getFileType
#
# @title "Gets the file type of a file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character in lower case letters.
# }
#
# \details{
#   By default, this methods returns the filename extension, but subclasses
#   may override this.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFileType", "AffymetrixFile", function(this, ...) {
  pattern <- "(.*)[.]([a-zA-Z0-9][a-zA-Z0-9]*)$";
  ext <- gsub(pattern, "\\2", this$.pathname);
  tolower(ext);
});


setMethodS3("getFileSize", "AffymetrixFile", function(this, ...) {
  file.info(this$.pathname)$size;
});


setMethodS3("fromFile", "AffymetrixFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (.checkArgs) {
    # Argument 'filename' and 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
  } else {
    pathname <- filename;
  }


  # Get all known subclasses (bottom up)
  clazz <- Class$forName(class(static)[1]);
  knownSubclasses <- rev(getKnownSubclasses(clazz));
  for (className in knownSubclasses) {
    clazz <- Class$forName(className);

    # Try reading the file using the static fromFile() method of each class
    static <- getStaticInstance(clazz);
    tryCatch({
      res <- fromFile(static, filename=pathname, .checkArgs=FALSE);
      return(res);
    }, error = function(ex) {})
  }

  throw("Could not read Affymetrix file '", pathname, "'.", 
             "Tried using all known subclasses of ", class(static)[1], ": ", 
                                     paste(knownSubclasses, collapse=", "));
}, static=TRUE)


setMethodS3("copyFile", "AffymetrixFile", function(this, filename, path=NULL, overwrite=FALSE, ..., verbose=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=!overwrite);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Assert that we're not trying to copy to itself
  if (identical(pathname, getPathname(this)))
    throw("Cannot copy Affymetrix file. Source and destination are identical: ", pathname);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Copy
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Copying Affymetrix file");
  verbose && cat(verbose, "Source: ", getPathname(this));
  verbose && cat(verbose, "Destination: ", pathname);
  # Get the pathname for the CEL file
  file.copy(getPathname(this), pathname, overwrite=overwrite);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create object of the same class.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- newInstance(this, pathname);

  res;
})


setMethodS3("stextSize", "AffymetrixFile", function(this, side=1, fmtstr="n=%d", size, pos=1, cex=0.7, col="darkgray", ...) {
  stext(side=side, line=-1, text=sprintf(fmtstr, round(size)), pos=pos, cex=cex, col=col, ...);
})

setMethodS3("stextLabel", "AffymetrixFile", function(this, side=3, fmtstr="%s", label=getLable(this), pos=0, cex=0.7, col="black", ...) {
  stext(side=side, text=sprintf(fmtstr, label), pos=pos, cex=cex, col=col, ...);
})

setMethodS3("stextLabels", "AffymetrixFile", function(this, others=NULL, side=3, fmtstr="%d) %s. ", pos=0, cex=0.7, col="black", ...) {
  # Build list of AffymetrixFile objects
  if (!is.list(others))
    others <- list(others);
  objects <- c(list(this), others);

  # Build text labels
  text <- vector("list", length(objects));
  for (kk in seq(along=objects)) {
    object <- objects[[kk]];
    if (is.null(object))
      next;
    value <- getLabel(object);
    str <- NULL;
    tryCatch({ str <- sprintf(fmtstr, value) }, error=function(ex){});
    if (is.null(str))
      str <- sprintf(fmtstr, kk, value);
    text[[kk]] <- str;
  }

  # Combine the into one
  text <- unlist(text, use.names=FALSE);
  text <- paste(text, collapse="");

  # Display it
  stext(side=side, text=text, pos=pos, cex=cex, col=col, ...);
})


############################################################################
# HISTORY:
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
