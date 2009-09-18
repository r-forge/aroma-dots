###########################################################################/**
# @RdocClass BaseFile
#
# @title "The BaseFile class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{sections}{A @list BaseFileSection objects.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
# @examples "../incl/BaseFile.Rex"
#
# @author
#*/###########################################################################
setConstructorS3("BaseFile", function(sections=NULL, ...) {
  if (is.list(sections)) {
    for (kk in seq(along=sections)) {
      # First, convert section to a BaseFileSection.
      section <- BaseFileSection(as.list(sections[[kk]]));

      # Then, try to convert it to a more specific subclass...
      type <- getType(section);
      type <- strsplit(type, split="[ ]+")[[1]];
      type <- capitalize(type);
      type <- paste(type, collapse="");
      className <- paste("BaseFile", type, sep="");
      clazz <- NULL;
      tryCatch({
        clazz <- Class$forName(className);
      }, error = function(ex) {
#         warning("BaseFileSection of type '", getType(section), 
#                 "' detected, but the corresponding class '", className, 
#                 "' was not found/implemented yet. Will use the superclass",
#                 " BaseFileSection instead.");
      })
      if (!is.null(clazz)) {
        section <- newInstance(clazz, as.list(section));
      }
      sections[[kk]] <- section;
    }
  }

  extend(Object(), "BaseFile",
    sections = sections
  )
})


#########################################################################/**
# @RdocMethod as.character
#
# @title "Gets a string description of object"
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
#   Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("as.character", "BaseFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- paste(class(this)[1], ":", sep="");
  s <- paste(s, " Number of sections: ", nbrOfSections(this), sep="");
  s;
})


#########################################################################/**
# @RdocMethod as.list
#
# @title "Gets a list representation of object"
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
#   Returns a @list.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("as.list", "BaseFile", function(x, ...) {
  # To please R CMD check...
  this <- x;

  lapply(this$sections, FUN=as.list);
})



#########################################################################/**
# @RdocMethod equals
#
# @title "Checks if this object equals another"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{other}{Another object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE if the objects are equal, otherwise @FALSE.
# }
#
# \details{
#  If the other object is not a BaseFile object, this method returns @FALSE.
#  Otherwise, \code{as.list()} is called on both objects are these are
#  compared with @see "base::identical".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("equals", "BaseFile", function(this, other, ...) {
  if (!inherits(other, "BaseFile"))
    return(FALSE);

  identical(as.list(this), as.list(other));
})



#########################################################################/**
# @RdocMethod nbrOfSections
# @aliasmethod length 
#
# @title "Gets the number of sections in this BASE file structure"
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
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("nbrOfSections", "BaseFile", function(this, ...) {
  length(this$sections);
})

setMethodS3("length", "BaseFile", function(x) {
  # To please R CMD check
  this <- x;

  length(this$sections);
}, appendVarArgs=FALSE)



#########################################################################/**
# @RdocMethod isSerial
#
# @title "Checks if a BASE file is in serial format or not"
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
#   Returns @TRUE if in serial format, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("isSerial", "BaseFile", function(this, ...) {
  spots <- getSections(this, "spots");
  assays <- lapply(spots, FUN=getAssays);
  counts <- lapply(assays, FUN=length);
  counts <- unlist(counts);

  isSerial <- all(counts <= 1);
  isSerial;
})


#########################################################################/**
# @RdocMethod getAllDataFiles
#
# @title "Gets the names of a external data files used by a BASE file"
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
#   Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getAllDataFiles", "BaseFile", function(this, ...) {
  # Currently, it is only 'spots' sections that "cache" data on file.
  spots <- getSections(this, "spots");
  files <- unlist(lapply(spots, FUN=getDataFiles));
  files <- unique(files);
  files;
})

#########################################################################/**
# @RdocMethod removeAllDataFiles
#
# @title "Removes all external data files used by a BASE file structure"
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
#   Returns (invisibly) @character @vector of the (successfully) 
#   removed files.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("removeAllDataFiles", "BaseFile", function(this, ...) {
  files <- getAllDataFiles(this);
  if (length(files) == 0)
    return(NULL);

  removedFiles <- c();
  for (file in files) {
    if (isFile(file)) {
      tryCatch({
        file.remove(file);
        removedFiles <- c(removedFiles, file);
      }, error = function(ex) {
      })
    }
  }

  failedFiles <- setdiff(files, removedFiles);
  if (length(failedFiles) > 0) {
    throw("Could not remove the following files: ", 
                                         paste(failedFiles, collapse=", "));
  }

  invisible(removedFiles);
})


#########################################################################/**
# @RdocMethod seq
#
# @title "Gets an index sequence for the sections in this BASE file structure"
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
#   Returns an @integer @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "nbrOfSections".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("seq", "BaseFile", function(this, ...) {
  seq(length=nbrOfSections(this));
})



#########################################################################/**
# @RdocMethod getSections
#
# @title "Gets a subset or all BASE sections"
# 
# \description{
#   @get "title" specified by their names and/or their header contents.
# }
#
# @synopsis
#
# \arguments{
#  \item{types}{A @vector of @character strings specifying which sections
#    to retrieve.}
#  \item{headers}{A @vector of @character strings specifying headers that
#    the sections must contain in order to be returned.}
#  \item{regexpr}{If @TRUE, the \code{types} and \code{headers} are 
#    interpreted as regular expressions, other exact matching is required.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @list of @see "BaseFileSection"s.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getSections", "BaseFile", function(this, types=NULL, headerValues=NULL, regexpr=FALSE, ...) {
  if (is.null(types)) {
    sections <- this$sections;
  } else {
    if (length(this$sections) == 0)
      return(NULL);
    sectionTypes <- names(this$sections);
    pattern <- types;
    if (!regexpr)
      pattern <- paste("^", pattern, "$", sep="");
    idx <- which(regexpr(pattern, sectionTypes) != -1);
    sections <- this$sections[idx];
  }

  if (length(sections) == 0)
    return(NULL);

  if (length(headerValues) > 0) {
    headerNames <- names(headerValues);
    keep <- c();
    for (kk in seq(along=sections)) {
      section <- sections[[kk]];
      ok <- TRUE;
      for (ll in seq(length(headerValues))) {
        headerName <- headerNames[ll];
        headerValue <- headerValues[[ll]];
        value <- section$headers[[headerName]];
        if (regexpr && !is.null(value) && !is.null(headerValue)) {
          ok <- ok && (regexpr(headerValue, value) != -1);
        } else {
          ok <- ok && (value %in% headerValue);
        }
        if (!ok)
          break;
      }
      if (ok)
        keep <- c(keep, kk);
    }
    sections <- sections[keep];
    if (length(sections) == 0)
      return(NULL);
  }

  sections;
})


#########################################################################/**
# @RdocMethod getSection
#
# @title "Gets one section of this BASE file structure"
# 
# \description{
#   @get "title".  If more than one exists of the same type, only the first
#   is returned.
# }
#
# @synopsis
#
# \arguments{
#  \item{type}{A @character string specifying the section to retrieve.}
#  \item{...}{Arguments passed to @seemethod "getSections".}
# }
#
# \value{
#   Returns a @see "BaseFileSection".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getSection", "BaseFile", function(this, type, ...) {
  type <- Arguments$getCharacter(type, nchar=c(1,Inf));
  getSections(this, types=type, ...)[[1]];
})



#########################################################################/**
# @RdocMethod appendSection
#
# @title "Adds a section to this BASE file structure"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{section}{A @see "BaseFileSection" to be appended.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE if section was added.
# }
#
# @author
#
# \seealso{
#   @seemethod "getSection".
#   @seemethod "removeSection".
#   @seemethod "replaceSection".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("appendSection", "BaseFile", function(this, section, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'section':
  if (!inherits(section, "BaseFileSection")) {
    throw("Argument 'section' is not a BaseFileSection object: ", 
                                                         class(section)[1]);
  }

  names <- c(names(this$sections), getType(section));
  this$sections <- append(this$sections, list(section));
  names(this$sections) <- names;

  invisible(TRUE);
})



#########################################################################/**
# @RdocMethod replaceSection
#
# @title "Replaces a section in a BASE file structure"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{oldSection}{The @see "BaseFileSection" to be removed.}
#  \item{newSection}{The @see "BaseFileSection" to replace the old one.}
#  \item{append}{A @logical. If \code{oldSection} is not found, should 
#     \code{newSection} be appended anyway?} 
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) @TRUE if section was replaced, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("replaceSection", "BaseFile", function(this, oldSection, newSection, append=FALSE, ...) {
  # Argument 'oldSection':
  if (!inherits(oldSection, "BaseFileSection"))
    throw("Argument 'oldSection' is not a BaseFileSection: ", class(oldSection)[1]);

  # Argument 'newSection':
  if (!inherits(newSection, "BaseFileSection"))
    throw("Argument 'newSection' is not a BaseFileSection: ", class(newSection)[1]);

  # Argument 'append':
  append <- Arguments$getLogical(append);
 

  sections <- this$sections;

  # Find position to replace
  replacePos <- NA;
  for (kk in seq(length=length(sections))) {
    if (equals(oldSection, sections[[kk]]))
      replacePos <- kk;
  }

  # If not found, what to do?
  if (is.na(replacePos)) {
    if (!append)
      return(invisible(FALSE));
    return(appendSection(this, newSection));
  }

  sections[[replacePos]] <- newSection;

  this$sections <- sections;

  invisible(TRUE);
})


#########################################################################/**
# @RdocMethod removeSection
#
# @title "Removes a section from a BASE file structure"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{section}{A @see "BaseFileSection" to be removed.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) the number of sections removed.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("removeSection", "BaseFile", function(this, section, ...) {
  # Argument 'section':
  if (!inherits(section, "BaseFileSection"))
    throw("Argument 'section' is not a BaseFileSection: ", class(section)[1]);


  sections <- this$sections;
  toRemove <- c();
  for (kk in seq(length=length(sections))) {
    if (equals(section, sections[[kk]]))
      toRemove <- c(toRemove, kk);
  }

  sections <- sections[-toRemove];
  this$sections <- sections;

  invisible(length(toRemove));
})



#########################################################################/**
# @RdocMethod read
#
# @title "Static method to read a BASE file structure"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{file}{A filename or a @connection to be read.}
#  \item{path}{If \code{file} is not a @connection, an optional path to
#    be added to the file to be read.}
#  \item{...}{Additional arguments passed to @see "readBaseFile", 
#    especially \code{extractSpotsData}.}
# }
#
# \value{
#   Returns an @see "BaseFile" object.
# }
#
# @author
#
# \seealso{
#   @seemethod "write".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("read", "BaseFile", function(static, file, path=NULL, ...) {
  if (!inherits(file, "connection")) {
    file <- Arguments$getReadablePathname(file, path, absolutePath=FALSE);

    if (regexpr("[.]gz$", file) != -1) {
      verbose && cat(verbose, "Extracting file: ", file);
      tmpname <- tempfile()
      n <- gunzip(file, tmpname);
      file <- tmpname;
      on.exit(file.remove(tmpname));
    }
  }

  # Read BASE file (or connection)
  sections <- readBaseFile(file, ...);

  # Return BaseFile object
  BaseFile(sections);
}, static=TRUE)



#########################################################################/**
# @RdocMethod write
#
# @title "Writes this BASE file structure"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{file}{A filename or a @connection to be written.}
#  \item{...}{If \code{file} is not a @connection, arguments passed to
#    to \code{\link[R.utils]{getWritablePathname.Arguments}()}, 
#    e.g. \code{path}, \code{overwrite}, and \code{mkdirs}.}
# }
#
# \value{
#   Returns what @see "writeBaseFile" returns.
# }
#
# @author
#
# \seealso{
#   @seemethod "read".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("write", "BaseFile", function(this, file, ...) {
  if (!inherits(file, "connection")) {
    file <- Arguments$getWritablePathname(file, ..., absolutePath=FALSE);
  }

  # Write BASE file
  writeBaseFile(con=file, as.list(this), ...);
})




#########################################################################/**
# @RdocMethod hasSerialBioAssaySets
#
# @title "Checks if the bioassay sets have data in serial format"
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
#   Returns a @logical value.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("hasSerialBioAssaySets", "BaseFile", function(this, ...) {
  spots <- getSections(this, "spots");
  nbrOfSections <- length(spots);
  if (nbrOfSections == 0)
    return(TRUE);
  for (kk in 1:nbrOfSections) {
    spot <- spots[[kk]];

    if (nbrOfAssays(spot) > 1)
      return(FALSE);
  }

  TRUE;
})


############################################################################
# HISTORY: 
# 2009-06-05
# o BUG FIX: getAllDataFiles() for BaseFile referenced object 'base' 
#   instead of 'this'.  Caught by the R code tools.
# 2005-12-21
# o Renamed argument 'name[s]' in getSection[s]() to 'type[s]' to be more
#   consistent with the 'type' argument of the BaseFileSection constructor.
# o BUG FIX: getSection() gave an "regexpr... invalid argument" error if
#   there were no sections.
# o BUG FIX: Yesterdays appendSection() did not work correctly and would
#   only add NA instead of the Object.
# 2005-12-20
# o Added appendSection().
# o replaceSection(old, new, append=TRUE), would not work if 'old' did not
#   exist.  Forgot to set name of new sections to its type.
# 2005-12-12
# o Updated the Rdoc comments for BaseFile$read().
# 2005-10-21
# o Updated the Rdoc comments to mention argument 'mustExist' and not
#   'mustExists'.
# 2005-07-20
# o BUG FIX: isSerial() used non-existing 'base' instead of 'this'.
# 2005-07-08
# o Added isSerial().
# 2005-07-07
# o Added replaceSection().
# 2005-06-19
# o Added Rdoc comments.
# o Added hasSerialBioAssaySets.
# 2005-06-07
# o Added Rdoc for equals().
# 2005-06-03
# o The constructor now tries to find the most specific BaseFileSection
#   class for each section. If no such is available, the BaseFileSection 
#   class is used.
# 2005-05-31
# o Added Rdoc comments.
# 2005-05-25
# o Created.
############################################################################
