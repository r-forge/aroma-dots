###########################################################################/**
# @RdocClass BaseFileSection
#
# @title "The BaseFileSection class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{section}{A @list structure.}
#   \item{type}{A @character string specifying the type of section.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
# 
# \details{
# }
#
#
# @author
#
# \seealso{
#   @see "BaseFile".
# }
#*/###########################################################################
setConstructorS3("BaseFileSection", function(section=NULL, type=NULL, ...) {
  if (inherits(section, "BaseFileSection"))
    return(clone(section));

  if (is.list(section)) {
    if (is.null(type)) {
      # Identify section type from headers
      type <- section$headers$section;
      if (is.null(type)) {
        throw("Argument 'type' is missing and argument 'section' contains no \"section\" header: ", paste(capture.output(str(section)), collapse="\n"));
      }
    }
  }

  section <- as.list(section);

  # Assert that there is a 'headers' structure.
  if (is.null(section$headers))
    section[["headers"]] <- list();

  # Update header 'section' according to type.
  if (is.null(section$headers$section))
    section$headers[["section"]] <- type;

  extend(Object(), "BaseFileSection",
    section = section,
    type = as.character(type)
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
setMethodS3("as.character", "BaseFileSection", function(this, ...) {
  s <- paste(class(this)[1], ":", sep="");
  s <- paste(s, " Type: ", getType(this), ".", sep="");
  s <- paste(s, " Number of headers: ", nbrOfHeaders(this), ".", sep="");
  if (hasData(this)) 
    s <- paste(s, " A data table exists.");
  s;
})



#########################################################################/**
# @RdocMethod equals
#
# @title "Checks if a BASE file section equals another"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{other}{Other object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @TRUE, if they are equal, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("equals", "BaseFileSection", function(this, other, ...) {
  if (!identical(class(this), class(other)))
    return(FALSE);

  if (!identical(getHeaders(this), getHeaders(other)))
    return(FALSE);

  if (!identical(getData(this), getData(other)))
    return(FALSE);

  TRUE;
})



#########################################################################/**
# @RdocMethod print
#
# @title "Prints a BASE file section to stdout"
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
#   Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("print", "BaseFileSection", function(x, ...) {
  # To please R CMD check
  this <- x;

  cat(as.character(this), "\n", sep="");
})



#########################################################################/**
# @RdocMethod validate
#
# @title "Validates BASE file section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{action}{A @character string specifying action if a problem is 
#               detected. If \code{"error"}, an exception is thrown,
#               if \code{"warning"}, a warning is generated, and if
#               \code{"none"}, to action is taken.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string explaining problems detected. If no problems
#   where detected, @NULL is returned.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("validate", "BaseFileSection", function(this, action=c("error", "warning", "none"), ...) {
  action <- match.arg(action);

  msg <- "";
  if (!hasData(this)) {
    if (length(getColumns(this)) != 0) {
      msg <- paste(msg, "No data table exists, but column names are specified: ", paste(getColumns(this), collapse=", "), ".", sep="");
    }
  }

  if (hasData(this)) {
    if (length(getColumns(this)) != sizeOfData(this)[2]) {
      msg <- paste(msg, " Number of columns names and number of columns in data table does not match: ", length(getColumns(this)), " != ", sizeOfData(this)[2], ".", sep="");
    }
  }

  msg <- trim(msg);
  if (nchar(msg) == 0)
    return(NULL);

  switch(action,
    error = throw(msg),
    warning = warning(msg)
  )

  msg;
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
setMethodS3("as.list", "BaseFileSection", function(x, ...) {
  # To please R CMD check...
  this <- x;

  this$section;
})


#########################################################################/**
# @RdocMethod getType
#
# @title "Gets the type of a BASE section"
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
setMethodS3("getType", "BaseFileSection", function(this, ...) {
  this$type;
})



#########################################################################/**
# @RdocMethod setType
#
# @title "Sets the type of a BASE section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{type}{A @character string specifying the new type.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) the old type as a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("setType", "BaseFileSection", function(this, type, ...) {
  oldType <- this$type;

  type <- Arguments$getCharacter(type, trim=TRUE, nchar=c(1,128));

  this$type <- type;
  invisible(oldType);
})




#########################################################################/**
# @RdocMethod nbrOfHeaders
#
# @title "Gets the number of headers in a BASE section"
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
setMethodS3("nbrOfHeaders", "BaseFileSection", function(this, ...) {
  length(this$section$headers);
})




#########################################################################/**
# @RdocMethod setHeader
#
# @title "Sets a header of a BASE section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{name}{A @character string of the name of the header.}
#  \item{value}{The value of the header. May be a @vector.}
#  \item{type}{A @character string specifying the new type.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) the old header value, or @NULL if missing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("setHeader", "BaseFileSection", function(this, name, value, ...) {
  name <- Arguments$getCharacter(name);

  # Assure that headers exists
  headers <- this$section$headers;
  if (is.null(headers))
    headers <- list();

  # Get old value
  oldValue <- headers[[name]];

  # Set new value
  headers[[name]] <- value;
  this$section$headers <- headers;

  invisible(oldValue);
})




#########################################################################/**
# @RdocMethod getHeaders
#
# @title "Gets the headers of a BASE section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{names}{A @vector of @character strings specifying which headers
#    to retrieve.}
#  \item{regexpr}{If @TRUE, the \code{names} are interpreted as regular 
#    expressions, other exact matching is required.}
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
#   @seemethod "getHeader".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getHeaders", "BaseFileSection", function(this, names=NULL, regexpr=FALSE, ...) {
  if (is.null(names))
    return(this$section$headers);

  headerNames <- names(this$section$headers);
  patterns <- names;
  if (!regexpr)
    patterns <- paste("^", patterns, "$", sep="");

  res <- list();

  for (kk in seq(along=patterns)) {
    pattern <- patterns[kk];
    idx <- which(regexpr(pattern, headerNames) != -1);
    value <- this$section$headers[idx];
    if (length(value) == 0) {
      value <- list(NULL);
      names(value) <- names[kk];
    }
    res <- c(res, value);
  }

  res;
})


#########################################################################/**
# @RdocMethod getHeader
#
# @title "Gets one header of a BASE section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{name}{A @character string specifying the header to retrieve.}
#  \item{...}{Arguments passed to @seemethod "getHeaders".}
# }
#
# \value{
#   Returns the value of a header.
# }
#
# @author
#
# \seealso{
#   @seemethod "getHeaders".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getHeader", "BaseFileSection", function(this, name, ...) {
  name <- Arguments$getCharacter(name);
  getHeaders(this, names=name, ...)[[1]];
})



#########################################################################/**
# @RdocMethod hasHeaders
#
# @title "Checks if specified headers exists in a BASE section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getHeaders".}
# }
#
# \value{
#   Returns a @logical @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "getHeaders".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("hasHeaders", "BaseFileSection", function(this, ...) {
  headers <- getHeaders(this, ...);
  if (length(headers) == 0)
    return(FALSE);

  isNull <- unlist(lapply(headers, FUN=is.null));
  !isNull;
})



#########################################################################/**
# @RdocMethod getColumns
#
# @title "Gets the column names of the data table"
# 
# \description{
#   @get "title", if it exists.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string @vector, if data table exists, otherwise
#   @NULL.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getColumns", "BaseFileSection", function(this, ...) {
  getHeader(this, "columns");
})


#########################################################################/**
# @RdocMethod hasColumns
#
# @title "Checks if certain columns exists or not"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{names}{A @character @vector of length one or more.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE for every column that exists, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("hasColumns", "BaseFileSection", function(this, names, ...) {
  names <- Arguments$getCharacters(names, length=c(1,Inf), nchar=c(1,Inf));
  (names %in% getColumns(this));
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
setMethodS3("setColumns", "BaseFileSection", function(this, names, ...) {
#  if (!hasData(this))
#    throw("Cannot set column names because there is no data table.");

  setHeader(this, "columns", names);
})



#########################################################################/**
# @RdocMethod addColumn
#
# @title "Adds a column to the other column names"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{name}{A column name.}
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
setMethodS3("addColumn", "BaseFileSection", function(this, name, ...) {
  name <- Arguments$getCharacter(name);

  columns <- getColumns(this);
  if (name %in% columns)
    throw("Argument 'name' is already an existing column: ", name);

  columns <- c(columns, name);
  setColumns(this, columns);
})



#########################################################################/**
# @RdocMethod hasData
#
# @title "Checks if section has a data table"
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
#   Returns @TRUE if a data table exists, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("hasData", "BaseFileSection", function(this, ...) {
  (length(this$section$data) > 0);
})



#########################################################################/**
# @RdocMethod setData
#
# @title "Sets data of a BASE file section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{data}{A @list structure.}
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
setMethodS3("setData", "BaseFileSection", function(this, data, ...) {
  # Argument 'data':
  if (is.matrix(data))
    data <- as.data.frame(data);

  if (!is.data.frame(data))
    throw("Argument 'data' is not a data.frame: ", mode(data));

  # Convert all factor fields into character fields.
  for (kk in seq(length=ncol(data))) {
    if (is.factor(data[[kk]]))
      data[[kk]] <- as.character(data[[kk]]);
  }

  this$section$data <- data;
})



#########################################################################/**
# @RdocMethod getData
#
# @title "Gets all or a subset of data of a BASE file section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{fields}{A @character @vector of field names to be retrieved.
#     If @NULL, all fields are returned.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @data.frame.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getData", "BaseFileSection", function(this, fields=NULL, ...) {
  getDataByKeys(this, byField=NULL, keys=NULL, fields=fields);
})


#########################################################################/**
# @RdocMethod getDataByKeys
#
# @title "Gets a subset of data for one or all bioassays in the bioassay set"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{byField}{A @character string of the field to be stratified on.
#     If @NULL, no stratification is done.}
#  \item{keys}{Rows with these values on the key field will be returned.}
#  \item{fields}{A @character @vector of the fields to be returned.
#     If @NULL, all fields are returned.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @data.frame.
# }
#
# @author
#
# \seealso{
#   @seemethod "getData".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getDataByKeys", "BaseFileSection", function(this, byField=NULL, keys=NULL, fields=NULL, ...) {
  # Argument 'byField':
  byField <- Arguments$getCharacter(byField, nchar=c(1,Inf));

  # Argument 'keys':
  keys <- Arguments$getVector(keys, length=c(0,Inf));

  # Argument 'field':
  fields <- Arguments$getCharacters(fields, length=c(0,Inf), nchar=c(1,Inf));

  # Nothing to do?
  if (!hasData(this))
    return(NULL);

  columns <- getColumns(this);
  if (!is.null(fields)) {
    fields <- as.character(fields);
    missing <- fields[!fields %in% columns];
    if (length(missing) > 0)
      throw("Unknown fields: ", paste(missing, collapse=", "));
  }

  data <- this$section$data;
  colnames(data) <- columns;

  if (is.null(byField)) {
    keep <- rep(TRUE, length=nrow(data));
  } else {
    if (!byField %in% columns)
      throw("Argument 'byField' is not an existing column: ", byField);

    keep <- rep(FALSE, length=nrow(data));
    for (kk in seq(length=length(keys))) {
      key <- keys[kk];
      keep <- keep | (data[[byField]] == key);
    }
  }

  if (is.null(fields)) {
    data <- data[keep,,drop=FALSE];
  } else {
    data <- data[keep,fields,drop=FALSE];
  }

  data;
})



#########################################################################/**
# @RdocMethod sizeOfData
#
# @title "Gets the dimension of the data table"
# 
# \description{
#   @get "title", if it exists.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   If data table exists, an @integer @vector of length two specifying the
#   number of rows and columns of the table is returned. Otherwise, @NULL 
#   is returned.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("sizeOfData", "BaseFileSection", function(this, ...) {
  if (!hasData(this))
    return(NULL);

  dim(getData(this));
})




#########################################################################/**
# @RdocMethod setDataFieldByKeys
#
# @title "Sets a data field by another key field"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{byField}{A string @character string of an existing column.}
#  \item{keys}{Rows with these values on the key field will be returned.}
#  \item{field}{A string @character string of an existing column.}
#  \item{...}{Arguments with names of keys and values to be assigned to the
#     field.}
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
setMethodS3("setDataFieldByKeys", "BaseFileSection", function(this, byField, keys, field, values, ...) {
  columns <- getColumns(this);

  # Argument 'byField':
  byField <- Arguments$getCharacter(byField);
  if (!byField %in% columns)
    throw("Argument 'byField' is not a existing column: ", byField);

  # Argument 'keys':
  keys <- Arguments$getVector(keys, length=c(1,Inf));

  # Argument 'field':
  field <- Arguments$getCharacter(field);

  # Argument 'values':
  values <- Arguments$getVector(values, length=c(1,1)*length(keys));

  # Assure that there is a <field> column, otherwise add one.
  if (!hasColumns(this, field))
    addColumn(this, field);

  # Nothing to do?
  if (length(keys) == 0)
    return();

  # Check for existance of key fields
  data <- this$section$data;
  missing <- keys[!keys %in% data[[byField]]];
  if (length(missing) > 0) {
    throw("Non-existing keys: ", paste(missing, collapse=", "), 
                       " not in ", paste(data[[byField]], collapse=", "));
  }

  for (kk in 1:length(keys)) {
    key <- keys[kk];
    value <- values[kk];
    if (length(value) != 1)
      throw("Argument 'values' should contain elements each of length one.");
    data[data[[byField]] == key, field] <- as.character(value);
  }

  # Update data table
  this$section$data <- data;
})


############################################################################
# HISTORY: 
# 2005-07-08
# o Added getDataByKeys().
# 2005-07-06
# o Now setHeader() makes sure 'headers' list structure exists.
# o Now setColumns() does not assert that a data table exists.
# o Now setData() assert that value is a data.frame or a matrix (not a 
#   list as before).
# o Added hasColumns(), setColumns() and addColumn().
# o Added setDataFieldByKeys().
# 2005-06-19
# o Added all Rdoc comments.
# 2005-06-15
# o Added hasData() and getColumns().
# 2005-06-07
# o Remove getData(). See BaseFileSpots class instead.
# 2005-05-31
# o Added Rdoc comments.
# 2005-05-25
# o Created.
############################################################################
