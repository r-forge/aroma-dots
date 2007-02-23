# Added '...' to some base functions. These will later be
# turned into default functions by setMethodS3().

colSums <- appendVarArgs(base::colSums)
colMeans <- appendVarArgs(base::colMeans)
write <- appendVarArgs(base::write)
append <- appendVarArgs(base::append)
getPackageName <- appendVarArgs(methods::getPackageName)


############################################################################
# HISTORY:
# 2007-02-23 [KS]
# o Make explicit reference to 'base' - this is safer, in case of colMeans()
#   defined by user or other packages.
# 2006-03-24
# o Created to please R CMD check.
############################################################################
