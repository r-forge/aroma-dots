
setMethodS3("getAbsolutePath", "MicroarrayData", function(this, filename, path=NULL) {
  warning("getAbsolutePath() in MicroarrayData is deprecated. Use package R.utils instead.");
  getAbsolutePath(filePath(path, filename));
}, protected=TRUE, static=TRUE, deprecated=TRUE);



setMethodS3("getField", "MicroarrayData", function(this, field, viewFunction=NULL) {
  # Assert that the specified 'field' is of type string.
  if (!is.character(field) || length(field) != 1)
    throw("Argument 'field' must be a single string.");

  # Assert that the specified field exists.
  if (!is.element(field, getFieldNames(this)))
    throw("No such field: ", field);

  attr(this, ".memberAccessorOrder") <- c(2,3,1,4,5);

  if (!is.null(viewFunction)) {
    if (!is.function(viewFunction))
      throw("Argument 'viewFunction' is not a function.");
    viewFunction(this[[field]])
  } else {
    this[[field]]
  }
}, private=TRUE, deprecated=TRUE)
    


#########################################################################/**
# @set "class=MicroarrayData"
# @RdocMethod set
#
# @title "Sets a subset of data specified by fields, slides and/or spot indices"
#
# @synopsis
#
# \arguments{
#   \item{value}{The values that should replace the specified subset of the
#     data. This argument should be a vector, matrix or a data frame and it
#     must have the same number of columns as the number of specified fields
#     (if fields are given).}
#   \item{fields}{The field names to be set. If \@NULL, all fields are set.}
#   \item{slides}{The slides that should be included. If \@NULL,
#     all slides are included.}
#   \item{index}{The spot indices that should be included. If \@NULL, 
#     all spots are included.}
#   Any of the above arguments are optional.
# }
#
# \description{
#   @get "title".
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   # Set the log ratios (M) in slide #3, to 1.0 for spot #5 and
#   # to 2.0 for spot #6:
#   set(ma, "M", slide=3, index=5:6, value=c(1.0,2.0))
#
#   # Get the log ratios (M) for spots 4-7 in the slide 2,3 and 4.
#   extract(ma, c("slide", "M"), slides=2:4, index=4:7)
#   # Gives:
#   #        slide           M
#   #  6388      2  0.78245098
#   #  6389      2  0.77639293
#   #  6390      2  0.93738005
#   #  6391      2 -0.07232743
#   #  12772     3  0.66280278
#   #  12773     3  1.00000000  <--
#   #  12774     3  2.00000000  <--
#   #  12775     3 -0.32318001
#   #  19156     4  1.21995595
#   #  19157     4  1.18084938
#   #  19158     4  1.32323496
#   #  19159     4 -0.04182559
# }
#
# @author
#
# \seealso{
#   @seemethod "extract".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("set", "MicroarrayData", function(this, value, fields=NULL, slides=NULL, index=NULL) {
  # Can't set the "slide" field. It is just a virtual field.
  fields <- setdiff(fields, "slide");

  # Gets the data frame with the specified fields.
  df <- extract(this, fields=c("slide", fields));

  # Get the rows that match the slides argument.
  if (!is.null(slides)) {
    slide.rows <- c();
    for (slide in slides) {
      slide.rows <- c(slide.rows, which(df$slide == slide));
    }
  } else {
    slide.rows <- NULL;
  }

  # Get the rows that match the index argument.
  if (!is.null(index)) {
    index.offsets <- (df$slide-1)*nbrOfSpots(this);
    index.rows <- unique(index.offsets + index);
  } else {
    index.rows <- NULL;
  }
  
  rows <- NULL;
  if (!is.null(slide.rows)) {
    if (!is.null(index.rows))
      rows <- intersect(slide.rows, index.rows)
    else
      rows <- slide.rows;
  } else if (!is.null(index.rows)) {
    rows <- index.rows;
  }

  if (is.null(rows) && is.null(fields)) {
    sub.df <- value;
  } else if (is.null(rows)) {
    df[,fields] <- value;
    sub.df <- df[,fields];
  } else {
    df[rows,fields] <- value;
    sub.df <- df[,fields];
  }

  # Now set all fields at once call the abstract method setField
  setField(this, value=sub.df, fields=fields);
  invisible(this);
}, private=TRUE, deprecated=TRUE);


############################################################################
# HISTORY:
# 2005-07-21
# o Made getAbsolutePath() of MicroarrayData deprecated. Alsom it is now 
#   internally using the same method in R.utils.
# 2005-05-03
# o Removed deprecated setField() since it was abstract and caused 
#   subclasses to be abstract, if not overriden.
# 2002-11-12
# o Created. Here all deprecated functions will be put in quarantine until
#   it is safe to remove them. 
############################################################################
