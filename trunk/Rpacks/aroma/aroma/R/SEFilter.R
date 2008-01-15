#########################################################################/**
# @RdocClass SEFilter
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
setConstructorS3("SEFilter", function(input, ...) {
  ok <- !missing(input);
  if (ok && !inherits(input, "MicroarrayData"))
    throw("Argument 'input' is not of class MicroarrayData: ", data.class(input));

  extend(FieldFilter(input, "SE", ...), "SEFilter"
  )
})

############################################################################
# HISTORY:
# 2003-04-21
# o Added Rdocs.
# 2002-02-26
# * Modified code to make use of setMethodS3's.
# 2001-09-29
# * Created by extracting the code from FieldFilter. Each class should be
#   in a seperate file!
############################################################################
