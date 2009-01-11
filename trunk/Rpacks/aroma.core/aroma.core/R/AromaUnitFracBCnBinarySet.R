###########################################################################/**
# @RdocClass AromaUnitFracBCnBinarySet
#
# @title "The AromaUnitFracBCnBinarySet class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitFracBCnBinarySet object represents a set of 
#  @see "AromaUnitFracBCnBinaryFile"s with \emph{identical} chip types.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaUnitSignalBinarySet".}
# }
#
# \details{
#   The term "allele B fraction" is also know as "allele B frequency", which
#   was coined by Peiffer et al. (2006).  Note that the term "frequency" is
#   a bit misleading since it is not a frequency in neither the statistical
#   nor the population sense, but rather only proportion relative to the
#   total amount of allele A and allele B signals, which is calculated for
#   each sample independently.
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \references{
#   \item Peiffer et al., \emph{High-resolution genomic profiling of 
#      chromosomal aberrations using Infinium whole-genome genotyping},
#      Genome Res, 2006.
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AromaUnitFracBCnBinarySet", function(...) {
  extend(AromaUnitSignalBinarySet(...), "AromaUnitFracBCnBinarySet");
})


setMethodS3("byName", "AromaUnitFracBCnBinarySet", function(static, name, tags=NULL, ..., chipType=NULL, paths="cnData") {
  suppressWarnings({
    path <- findByName(static, name=name, tags=tags, chipType=chipType, 
                                           ..., paths=paths, mustExist=TRUE);
  })

  suppressWarnings({
    fromFiles(static, path=path, ..., pattern=".*,(frac|freq)B[.]asb$");
  })
}, static=TRUE) 




############################################################################
# HISTORY:
# 2009-01-03
# o Renamed from freqB to fracB, because it is a fraction, not a frequency.
# 2008-05-11
# o Created.
############################################################################
