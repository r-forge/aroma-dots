###########################################################################/**
# @RdocClass RmaCnPlm
#
# @title "The RmaCnPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the log-additive model used in RMA.
#  It can be used to fit the model on a @see "AffymetrixCelSet".
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RmaSnpPlm".}
#   \item{combineAlleles}{If @FALSE, allele A and allele B are treated 
#      seperately, otherwise together.}
#   \item{name}{The name of the PLM, also used as part of the pathname.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#   TO DO.
# }
#
# @author
#
#*/###########################################################################
setConstructorS3("RmaCnPlm", function(..., combineAlleles=FALSE, name="modelCnRmaPlm") {
  extend(RmaSnpPlm(..., name=name), c("RmaCnPlm", uses(CnPlm())),
    combineAlleles = combineAlleles
  )
})


############################################################################
# HISTORY:
# 2006-09-10
# o Recreated.
############################################################################
