# Added '...' to some base functions. These will later be
# turned into default functions by setMethodS3().

readLines <- appendVarArgs(readLines);


############################################################################
# HISTORY:
# o Created to please R CMD check.
############################################################################
