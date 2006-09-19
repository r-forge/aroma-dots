###########################################################################/**
# @RdocClass CnChipEffectFile
#
# @title "The CnChipEffectFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in a copy-number probe-level
#  models.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "SnpChipEffectFile".}
#   \item{combineAlleles}{A @logical indicating if the signals from allele A and
#      allele B are combined or not.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#   An object of this class is typically part of a @see "CnChipEffectSet".
# }
#*/###########################################################################
setConstructorS3("CnChipEffectFile", function(..., combineAlleles=FALSE) {
  extend(SnpChipEffectFile(...), "CnChipEffectFile",
    combineAlleles = combineAlleles
  )
})

setMethodS3("getCellIndices", "CnChipEffectFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- digest(list(units=units, ...));
  res <- this$.cellIndices[[key]];
  if (!is.null(res)) {
    verbose && cat(verbose, "getCellIndices.CnChipEffectFile(): Returning cached data");
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get and restructure cell indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cells <- NextMethod("getCellIndices", this, ..., verbose=verbose);

  # If combining alleles, return only every second group.
  # In order to improve readability we merge the names of alleles groups
  # combined, e.g. groups 'C' and 'G' become group 'CG'.
  if (this$combineAlleles) {
    cells <- applyCdfGroups(cells, function(groups) {
      ngroups <- length(groups);
      odds <- seq(from=1, to=ngroups, by=2);
      names <- names(groups);
      groups <- groups[odds];
      if (ngroups >= 2) {
        evens <- seq(from=2, to=ngroups, by=2);
        names <- paste(names[odds], names[evens], sep="");
        names(groups) <- names;
      }
      groups;
    })
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "getCellIndices.CnChipEffectFile(): Updating cache");
  this$.cellIndices <- list();
  this$.cellIndices[[key]] <- cells;

  cells;
})


############################################################################
# HISTORY:
# 2006-09-17
# o Added an in-memory cache for getCellIndices().
# 2006-09-12
# o Updated and probably working. When combining alleles, the names of the
#   groups returned consist of the allele A and allele group names.
# 2006-09-11
# o Created.
############################################################################
