###########################################################################/**
# @RdocClass BaseFileSpots
#
# @title "The BaseFileSpots class extending BaseFileSection"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to constructor of superclass.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
# @author
#*/###########################################################################
setConstructorS3("BaseFileSpots", function(...) {
  extend(BaseFileSection(..., type="spots"), "BaseFileSpots",
    assayFieldSep="_of_"
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
setMethodS3("as.character", "BaseFileSpots", function(this, ...) {
  s <- paste(class(this)[1], ":", sep="");

  s <- paste(s, " Type: ", getType(this), sep="");

  assays <- getAssays(this);
  s <- paste(s, ". Number of assays: ", length(assays), sep="");

  if (length(assays) <= 4) {
    assays <- paste(assays, collapse=", ");
  } else {
    assays <- c(assays[c(1,2)], "...", assays[length(assays)]);
  }
  s <- paste(s, " (", paste(assays, collapse=", "), ")", sep="");

  # Try to get the layout
  tryCatch({
    layout <- getLayout(this);
    s <- paste(s, ". ", as.character(layout), sep="");
  }, error = function(ex) {});

  s <- paste(s, ". Number of headers: ", nbrOfHeaders(this), sep="");

  s <- paste(s, ". Columns: ", 
                            paste(getColumns(this), collapse=", "), sep="");
  s <- paste(s, ". Assay fields: ", 
                        paste(getAssayFields(this), collapse=", "), sep="");
  s;
})



#########################################################################/**
# @RdocMethod getAssays
#
# @title "Gets the ids of the assays in this section"
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
#   Returns a @character string @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getAssays", "BaseFileSpots", function(this, ...) {
  getHeader(this, "assays");
})



#########################################################################/**
# @RdocMethod setAssays
#
# @title "Sets the ids of the assays in this section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{ids}{The assay ids.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("setAssays", "BaseFileSpots", function(this, ids, ...) {
  ids <- Arguments$getVector(ids, length=c(1,Inf));
  setHeader(this, "assays", ids);
})



#########################################################################/**
# @RdocMethod nbrOfAssays
#
# @title "Gets the number of assays in this section"
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
setMethodS3("nbrOfAssays", "BaseFileSpots", function(this, ...) {
  length(getAssays(this));
})



#########################################################################/**
# @RdocMethod getColumns
#
# @title "Gets the names of the non-assay-field data columns"
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
#   Returns a @character string @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "getAssayFields"
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getColumns", "BaseFileSpots", function(this, ...) {
  columns <- NextMethod("getColumns");
  setdiff(columns, "assayData");
})


#########################################################################/**
# @RdocMethod setColumns
#
# @title "Sets the column names of the data table"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{names}{A @character @vector.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("setColumns", "BaseFileSpots", function(this, names, ...) {
  names <- Argument$getCharacters(names);

  columns <- NextMethod("getColumns");
  hasAssayData <- any("assayData" %in% columns);
  if (hasAssayData)
    names <- c(names, "assayData");

  NextMethod("setColumns", this, names);  
})



#########################################################################/**
# @RdocMethod getAssayFields
#
# @title "Gets the names of the assay-field data columns"
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
#   Returns a @character string @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "getColumns"
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getAssayFields", "BaseFileSpots", function(this, ...) {
  getHeader(this, "assayFields");
})


#########################################################################/**
# @RdocMethod getDataFiles
#
# @title "Gets the filenames where data of an assay-field section is stored"
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
#   Returns a @character @vector, or @NULL, if no data is stored on file.
# }
#
# @author
#
# \seealso{
#   @seemethod "setDataFiles"
#   @seemethod "hasDataFiles"
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getDataFiles", "BaseFileSpots", function(this, ...) {
  getHeader(this, "dataFiles");
})



#########################################################################/**
# @RdocMethod setDataFiles
#
# @title "Sets the filenames where data of an assay-field section is stored"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{filenames}{A @character @vector of filenames.}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "getDataFiles"
#   @seemethod "hasDataFiles"
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("setDataFiles", "BaseFileSpots", function(this, filenames, ...) {
  if (!is.null(filenames))
    filenames <- as.character(filenames);
  setHeader(this, "dataFiles", filenames);
})


#########################################################################/**
# @RdocMethod hasDataFiles
#
# @title "Checks if data of an assay-field section is stored on file"
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
#   @seemethod "getDataFiles"
#   @seemethod "setDataFiles"
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("hasDataFiles", "BaseFileSpots", function(this, ...) {
  length(getDataFiles(this) > 0);
})




#########################################################################/**
# @RdocMethod getData
#
# @title "Gets a subset of or the complete data table"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{assays}{A @character string @vector of assays to be returned.}
#  \item{fields}{A @character string @vector of fields, both non-assay fields ("columns") and assay fields to be returned.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @data.frame.
# }
#
# \details{
#  If only one assay is returned, the assay name suffices are excluded 
#  from the assay field columns, otherwise they are included.
# }
#
# @author
#
# \seealso{
#   @seemethod "getColumns"
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getData", "BaseFileSpots", function(this, assays=NULL, fields=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'assays':
  assays <- Arguments$getCharacters(assays, nchar=c(1,128));
  if (is.null(assays)) {
    assays <- getAssays(this);
  } else {
    missing <- assays[!(assays %in% getAssays(this))];
    if (length(missing) > 0) {
      throw("Argument 'assays' contained unknown assays: ", 
                                        paste(missing, collapse=", "));
    }
  }

  # Argument 'fields':
  fields <- Arguments$getCharacters(fields, nchar=c(1,128));


  # Split 'fields' into "column" and "assay" fields.
  columnFields <- getColumns(this);
  assayFields <- getAssayFields(this);
  if (!is.null(fields)) {
    columnFields <- fields[fields %in% columnFields];
    assayFields <- fields[fields %in% assayFields];
    missing <- setdiff(fields, c(columnFields, assayFields));
    if (length(missing) > 0) {
      throw("Argument 'fields' contained unknown columns and assay fields: ", 
                                        paste(missing, collapse=", "));
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data files or is it already read?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (hasDataFiles(this)) {
    data <- NULL;
    for (dataFile in getDataFiles(this)) {
      # Read all data...
      tmp <- read.table(file=dataFile, header=TRUE, sep="\t");
      # ...but keep only wanted fields
      if (is.null(data)) {
        data <- tmp[, c(columnFields, assayFields)];
      } else {
        tmp <- tmp[, assayFields];
        data <- cbind(data, tmp);
      }
      rm(tmp);
    }

    if (length(assays) > 1) {
      fieldnames <- columnFields;
      for (assay in assays) {
        fieldnames <- c(fieldnames, paste(assayFields, this$assayFieldSep, assay, sep=""));
      }
      colnames(data) <- fieldnames;
    }
  } else {
    # Get data
    data <- this$section$data;
  
    # Column fields
    colnames <- columnFields;
    
    # Assay fields
    if (length(assayFields) > 0) {
      for (assay in assays) {
        # First, try to get the fields, with assay name suffices...
        fieldnames <- paste(assayFields, this$assayFieldSep, assay, sep="");
    
        # Second, try without assay name suffices...
        if (!all(fieldnames %in% names(data)))
          fieldnames <- assayFields;
    
        colnames <- c(colnames, fieldnames);
      }
    }
  
    missing <- colnames[!(colnames %in% colnames(data))];
    if (length(missing) > 0) {
      throw("Internal exception: Seems that the read data table does not contain all fields as expected: ", paste(missing, collapse=", "));
    }
  
    data <- data[,colnames,drop=FALSE];
  
    if (length(assays) == 1 && length(assayFields) > 0) {
      # Remove assay suffix.
      pattern <- paste(this$assayFieldSep, assays, sep="");
      colnames(data) <- gsub(pattern, "", colnames);
    }
  }

  data;
})



#########################################################################/**
# @RdocMethod setDataFields
#
# @title "Sets data field for an assay in a BASE file spots section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{assay}{A @character string for an existing assay.}
#  \item{values}{A (named) @data.frame or a @list containing values to be 
#                                 assigned to the fields in \code{fields}.}
#  \item{fields}{A @character of fields to be assigned.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) the new set of data as a @data.frame.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("setDataFields", "BaseFileSpots", function(this, assay=NULL, values, fields=names(values), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'assay':
  assay <- Arguments$getCharacter(assay, length=c(0,1), nchar=c(1,128));
  if (is.null(assay)) {
    if (nbrOfAssays(this) > 1)
      throw("Argument 'assay' must be given when there are more than one assay in the 'spots' section: ", nbrOfAssays(this));
    assay <- getAssays(this)[1];
  } else {
    if (!assay %in% getAssays(this))
      throw("Argument 'assay' is not a known assay: ", assay);
  }

  # Argument 'values':
  values <- as.data.frame(values);
  if (length(values) == 0)
    throw("Argument 'values' has no fields.");

  # Argument 'fields':
  fields <- Arguments$getCharacters(fields, length=length(values), nchar=c(1,128));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Gets and validate field names
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Split 'fields' into "column" and "assay" fields.
  columnFields <- getColumns(this);
  assayFields <- getAssayFields(this);
  if (!is.null(fields)) {    
    columnFields <- fields[fields %in% columnFields];
    assayFields <- fields[fields %in% assayFields];
    missing <- setdiff(fields, c(columnFields, assayFields));
    if (length(missing) > 0) {
      throw("Argument 'fields' contained unknown columns and assay fields: ", 
                                        paste(missing, collapse=", "));
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assign to existing data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (hasDataFiles(this)) {
    # Read data from file
    filename <- getDataFiles(this)[1];
    data <- read.table(file=filename, header=TRUE, sep="\t");

    # Update/add new fields
    data[,fields] <- values;

    # Save data back to file
    write.table(data, file=filename, quote=FALSE, sep="\t", 
                             row.names=FALSE, col.names=TRUE);
  } else {
    # Get data
    data <- this$section$data;
  
    # Column fields
    colnames <- columnFields;
    
    # Assay fields
    if (length(assayFields) > 0) {
      for (assay in assays) {
        # First, try to get the fields, with assay name suffices...
        fieldnames <- paste(assayFields, this$assayFieldSep, assay, sep="");
    
        # Second, try without assay name suffices...
        if (!all(fieldnames %in% names(data)))
          fieldnames <- assayFields;
    
        colnames <- c(colnames, fieldnames);
      }
    } 

    # Update/add new fields
    data[,colnames] <- values;

    # Store updated data
    this$section$data <- data;
  }

  invisible(data);
})


############################################################################
# HISTORY: 
# 2005-12-20
# o BUG FIX: When reading data cached to file, data cells containing spaces
#   would generate an error in the internal read.table(). Forgot to specify
#   sep="\t" when reading cached data.
# 2005-09-06
# o Columns and assay fields are now reported by as.character().
# 2005-07-06
# o BUG FIX: getData() and setDataField() generated data file names on the
#   fly, but should really use getDataFiles(), e.g. the 'assays' header may
#   have changed.
# o Added setColumns(), but it is not tested.
# 2005-06-19
# o Added all Rdoc comments.
# 2005-06-18
# o Added setDataFields().
# 2005-06-15
# o Now getData() reads extracted data files too.
# 2005-06-07
# o Added Rdoc comments.
# o Renamed getDataFields() to getData().
# 2005-06-03
# o Modified getLayout(), getSpotPosition(), getForeground() and 
#   getBackground() such that they can retrieve one or multiple assays.
# o Added getDataFields().
# o Added getAssays() and nbrOfAssays().
# 2005-06-02
# o Created.
############################################################################
