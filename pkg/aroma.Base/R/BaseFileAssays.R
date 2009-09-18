###########################################################################/**
# @RdocClass BaseFileAssays
#
# @title "The BaseFileAssays class extending BaseFileSection"
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
setConstructorS3("BaseFileAssays", function(...) {
  extend(BaseFileSection(..., type="assays"), "BaseFileAssays")
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
setMethodS3("as.character", "BaseFileAssays", function(x, ...) {
  # To please R CMD check
  this <- x;
  getLayout <- NULL; rm(getLayout);

  s <- paste(class(this)[1], ":", sep="");
  s <- paste(s, " Type: ", getType(this), sep="");
  assays <- getName(this);
  s <- paste(s, ". Number of assays: ", nrow(assays), sep="");
  if (length(assays) <= 4) {
    assays <- paste(assays, collapse=", ");
  } else {
    assays <- c(assays[c(1,2)], "...", assays[length(assays)]);
  }
  s <- paste(s, " (", paste(assays, collapse=", "), ")", sep="");
  tryCatch({
    layout <- getLayout(this);
    s <- paste(s, ". ", as.character(layout), sep="");
  }, error = function(ex) {});
  s <- paste(s, ". Number of headers: ", nbrOfHeaders(this), sep="");
  s;
})


#########################################################################/**
# @RdocMethod validate
#
# @title "Validates assays section"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{action}{A @character string specifying action if a problem is 
#               detected. See @see "validate.BaseFileSection" for details.}
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
setMethodS3("validate", "BaseFileAssays", function(this, action=c("error", "warning", "quiet"), ...) {
  action <- match.arg(action);

  msg <- NextMethod("validate", action="none");

  if (!hasData(this)) {
    if (nbrOfAssays(this) != 0) {
      msg <- paste(msg, " No data table exists, but expected ", nbrOfAssays(this), " columns.", sep="");
    }
  }

  if (hasData(this)) {
    if (nbrOfAssays(this) != sizeOfData(this)[1]) {
      msg <- paste(msg, "Number of expected assays and number of rows in data table does not match: ", nbrOfAssays(this), " != ", sizeOfData(this)[1], ".", sep="");
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
}, protected=TRUE)



#########################################################################/**
# @RdocMethod getCount
#
# @title "Gets the (expected) number of assays in the data table"
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
#   Returns an @integer, or @NULL if count is not specified.
# }
#
# @author
#
# \seealso{
#   @seemethod "setCount"
#   @seemethod "nbrOfAssays"
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("getCount", "BaseFileAssays", function(this, ...) {
  count <- getHeader(this, "count");
  if (length(count) == 0)
    return(NULL);
  count <- as.integer(count);
  count;
})


#########################################################################/**
# @RdocMethod setCount
#
# @title "Sets (or update) the number of assays in the data table"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{count}{An @integer or @NULL. If @NULL, the number of assays is
#    set to the number of rows in the current data table.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer for the old count, or @NULL if missing.
# }
#
# @author
#
# \seealso{
#   @seemethod "getCount"
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("setCount", "BaseFileAssays", function(this, count=NULL, ...) {
  oldCount <- getCount(this);

  # Infer 'count' from data table?
  if (is.null(count)) {
    if (!hasData(this))
      throw("Argument 'count' is NULL and no data table exists.");
    count <- nrow(getData(this));
  }

  # Set new count
  count <- Arguments$getInteger(count);
  setHeader(this, "count", count);

  invisible(oldCount);
})


#########################################################################/**
# @RdocMethod nbrOfAssays
#
# @title "Gets the (expected) number of assays in the data table"
# 
# \description{
#   @get "title".  This method return what @seemethod "getCount" returns,
#   except that if the latter returns @NULL, this returns zero.
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
#   @seemethod "getCount"
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("nbrOfAssays", "BaseFileAssays", function(this, ...) {
  count <- getCount(this);
  if (length(count) == 0)
    count <- as.integer(0);
  count;
})




#########################################################################/**
# @RdocMethod getAnnotationColumns
#
# @title "Gets the names of the annotation columns"
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
setMethodS3("getAnnotationColumns", "BaseFileAssays", function(this, ...) {
  getHeader(this, "annotationColumns");
})


#########################################################################/**
# @RdocMethod setAnnotationColumns
#
# @title "Sets the names of the annotation columns"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{names}{A @character @vector of column names. 
#    If there already are existing annotation columns specified, the
#    number of new column names must match the existing number, otherwise
#    an exception is thrown.}
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
#   @seemethod "getAnnotationColumns"
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("setAnnotationColumns", "BaseFileAssays", function(this, names, ...) {
  # Argument 'names':
  names <- Arguments$getCharacters(names, nchar=c(1,Inf));

  # If there are already annotation columns, make sure the new ones are
  # equally many.
  currNames <- getAnnotationColumns(this);
  if (length(currNames) > 0) {
    if (length(currNames) != length(names)) {
      throw("The number of columns in argument 'names' does not match the number of existing annotation columns: ", length(names), " != ", length(currNames));
    }
  }

  setHeader(this, "annotationColumns", names);
})



#########################################################################/**
# @RdocMethod setData
#
# @title "Sets data of a BASE file section"
# 
# \description{
#   @get "title". This calls the samemethod in the superclass, and updates
#   the headers 'columns' and 'count' accordingly.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Passed to 
#            \code{\link[aroma.Base:setData.BaseFileSection]{setData()}}.}
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
setMethodS3("setData", "BaseFileAssays", function(this, ...) {
  NextMethod("setData");
  # Register column names here, because these are needed in setCount().
  setColumns(this, colnames(this$section$data));
  setCount(this);
})




#########################################################################/**
# @RdocMethod getDataByIds
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
#  \item{ids}{Ids of assays to be returned.}
#  \item{...}{Arguments passed to 
#  \code{\link[aroma.Base:getDataByKeys.BaseFileSection]{getDataByKeys()}}.}
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
setMethodS3("getDataByIds", "BaseFileAssays", function(this, ids, ...) {
  if (length(ids) == 0)
    return(NULL);

  getDataByKeys(this, byField="id", keys=ids, ...);
})



#########################################################################/**
# @RdocMethod setDataFieldByIds
#
# @title "Sets a field of one or all bioassays in the bioassay set"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments with names of assay ids and values to be assigned
#     to the field.}
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
setMethodS3("setDataFieldByIds", "BaseFileAssays", function(this, ids, field, values, ...) {
  if (length(ids) != length(values)) {
    throw("Argument 'ids' and 'values' are of different lengths: ", 
                                       length(ids), " != ", length(values));
  }

  setDataFieldByKeys(this, byField="id", keys=ids, field=field, values=values, ...);
})




#########################################################################/**
# @RdocMethod getId
#
# @title "Gets the id of all bioassays in the bioassay set"
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
setMethodS3("getId", "BaseFileAssays", function(this, ...) {
  value <- getData(this, fields="id")[[1]];
  as.character(value);
})



#########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of all bioassays in the bioassay set"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{ids}{Ids of assays to be returned. If @NULL, all are returned.}
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
setMethodS3("getName", "BaseFileAssays", function(this, ids=NULL, ...) {
  if (is.null(ids)) {
    value <- getData(this, fields="name")[[1]];
  } else {
    value <- getDataByIds(this, ids=ids, fields="name")[[1]];
  }
  as.character(value);
})



#########################################################################/**
# @RdocMethod setName
#
# @title "Sets the names of one or all bioassays in the bioassay set"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments with names of assay ids and values that are 
#     single @character strings of the corresponding assay names.}
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
setMethodS3("setName", "BaseFileAssays", function(this, ...) {
  setDataFieldByIds(this, field="name", ...);
})




#########################################################################/**
# @RdocMethod getParents
#
# @title "Gets the ids of the assay parents of all bioassays in the bioassay set"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{ids}{Ids of assays to be returned. If @NULL, all are returned.}
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
setMethodS3("getParents", "BaseFileAssays", function(this, ids=NULL, ...) {
  columns <- getColumns(this);
  if (!"parent" %in% columns)
    return(NULL);
  if (is.null(ids)) {
    value <- getData(this, fields="parent")[[1]];
  } else {
    value <- getDataByIds(this, ids=ids, fields="parent")[[1]];
  }
  as.character(value);
})


#########################################################################/**
# @RdocMethod setParents
#
# @title "Sets the ids of the assay parents of one or all bioassays in the bioassay set"
# 
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments with names of assay ids and values that are 
#     @character strings of parents, either as a @character @vector or
#     a singe @character string where parents are split by \code{/} 
#     (or \code{,} or \code{;}).}
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
setMethodS3("setParents", "BaseFileAssays", function(this, ids, parents, ...) {
  # Be forgiving if 'ids' is a single id and 'parents' has several elements.
  if (length(ids) == 1 && length(parents) > 1)
    parents <- list(parents);

  # Argument 'parents':
  parents <- as.list(parents);

  # Split and rejoin values by '/'.
  for (kk in seq(length=length(parents))) {
    value <- as.character(parents[[kk]]);

    # Split all parents
    value <- unlist(strsplit(value, split="[/,;]"));

    # Assert that no duplicated parents exists
    if (any(duplicated(value))) {
      throw("Duplicated parents detected: ", 
                           paste(value[duplicated(value)], collapse=", "));
    }

    # Join to internal format
    value <- paste(value, collapse="/");

    # Update list
    parents[[kk]] <- value;
  }

  parents <- unlist(parents);

  setDataFieldByIds(this, ids=ids, field="parents", values=parents);
})






############################################################################
# HISTORY: 
# 2005-12-20
# o Added setAnnotationColumns().
# 2005-07-08
# o Added optional argument 'ids' to getName() and getParents().
# o Added setData(), which also updates 'count' and 'columns'.
# 2005-07-06
# o Added setParents(), and setName().
# o Added generic setDataFieldByIds().
# 2005-06-19
# o Added all Rdoc comments.
# 2005-06-15
# o Created.
############################################################################
