###########################################################################/**
# @RdocClass GenotypeCallXdrSet
#
# @title "The GenotypeCallXdrSet class"
#
# \description{
#  @classhierarchy
#
#  The GenotypeCallXdrSet class represents a GenotypeCallSet where the
#  data files are stored in the binary XDR file format.
#
#  \emph{You should consider this class temporary until another file format
#  is defined for genotype calls. /HB 2006-12-14.}
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "GenotypeCallSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
# @visibility "private"
#*/###########################################################################
setConstructorS3("GenotypeCallXdrSet", function(...) {
  extend(GenotypeCallSet(...), "GenotypeCallXdrSet")
})


setMethodS3("fromFiles", "GenotypeCallXdrSet", function(static, ..., 
                   pattern=",calls[.]xdr$", fileClass="GenotypeCallXdrFile") {
  # S3 method dispatch does not work for static methods.
  fromFiles.GenotypeCallSet(static, ..., pattern=pattern, fileClass=fileClass);
}, static=TRUE)



##############################################################################
# HISTORY:
# 2006-12-14
# o Created.
##############################################################################
