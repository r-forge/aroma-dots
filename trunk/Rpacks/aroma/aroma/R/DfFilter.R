#########################################################################/**
# @RdocClass DfFilter
#
# @title "FieldFilter for the log-intensities (A)"
#
# \description{
#  @classhierarchy
#
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{input}{The input @see "MicroarrayData" object.}
#   \item{...}{Any arguments accepted by the @see "FieldFilter" constructor.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
#
# @author
#
# \seealso{
#   See also the @see "FieldFilter" class.
# }
#
# \keyword{manip}
#*/#########################################################################
setConstructorS3("DfFilter", function(input, ...) {
  ok <- !missing(input);
  if (ok && !inherits(input, "MicroarrayData"))
    throw("Argument 'input' is not of class MicroarrayData: ", data.class(input));

  extend(FieldFilter(input, "df", ...), "DfFilter"
  )
})

############################################################################
# HISTORY:
# 2003-04-21
# o Added Rdocs.
# 2002-02-26
# * Rewritten for setClassS3.
# 2001-09-29
# * Created by extracting the code from FieldFilter. Each class should be
#   in a seperate file!
############################################################################
