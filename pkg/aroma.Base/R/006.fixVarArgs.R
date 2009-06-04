# Added '...' to some base functions. These will later be
# turned into default functions by setMethodS3().

write <- appendVarArgs(write)


############################################################################
# HISTORY:
# 2005-06-19
# o Created to please R CMD check.
############################################################################
