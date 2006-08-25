# Added '...' to some base functions. These will later be
# turned into default functions by setMethodS3().

colSums <- appendVarArgs(colSums)
colMeans <- appendVarArgs(colMeans)
write <- appendVarArgs(write)


############################################################################
# HISTORY:
# 2006-03-24
# o Created to please R CMD check.
############################################################################
