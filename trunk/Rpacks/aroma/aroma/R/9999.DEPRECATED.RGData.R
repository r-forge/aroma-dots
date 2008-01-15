

#########################################################################/**
# @set "class=RGData"
# @RdocMethod setField
#
# @title "Sets one or many fields"
#
# @synopsis
#
# \arguments{
#   \item{value}{The values that should replace the specified subset of the
#     data. This argument should be a vector, matrix or a data frame and it
#     must have the same number of columns as the number of specified fields.}
#   \item{fields}{The field names to be set. If \@NULL, all fields are
#     set. Valid values are \code{"G"} and \code{"R"}.}
# }
#
# \description{
#   Sets one or many fields. The field names must be specified by the 
#   argument \code{field} or if \code{field} is \@NULL all fields are
#   considered. The argument \code{value} must be able to be converted to
#   a matrix by \code{as.matrix(value)} where the resulting matrix must
#   have the same number of columns as the specified number of fields.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#   rg <- as.RGData(ma)
#
#   # Set all R's to zeros.
#   setField(rg, "R", rep(0, rg$size()))
#   range(extract(rg, "R"))  # 0 0
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("setField", "RGData", function(this, fields, value) {
  value <- as.matrix(value);
  ncol <- nbrOfSlides(this);

  # df now contains the final settings. Now we just update field by field
  # without having to think about the arguments slides and index.
  k <- 1;
  for (field in fields) {
    if (field == "R" || field == "G") {
      this[[field]] <- matrix(value[,k], ncol=ncol);
    } else {
      throw("Trying to set field \"", field, "\". Fields that can be set are \"R\" and \"G\".");
    }
    k <- k + 1;
  }
  invisible(this);
}, private=TRUE, deprecated=TRUE);
            
            


############################################################################
# HISTORY:
# 2002-11-12
# o Created. Here all deprecated functions will be put in quarantine until
#   it is safe to remove them. 
############################################################################
