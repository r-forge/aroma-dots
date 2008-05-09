# Added '...' to some base functions. These will later be
# turned into default functions by setMethodS3().

colSums <- appendVarArgs(colSums);
colMeans <- appendVarArgs(colMeans);
write <- appendVarArgs(write);
append <- appendVarArgs(append);
getPackageName <- appendVarArgs(getPackageName);
readLines <- appendVarArgs(readLines);


############################################################################
# HISTORY:
# 2007-02-27 [HB]
# o BUG FIX: Removed explicit reference to 'base' etc again. The reason is 
#   that if a previous package already modified, say, write(), to become a 
#   generic function, that was overwritten again when this package was 
#   loaded.
# 2007-02-23 [KS]
# o Make explicit reference to 'base' - this is safer, in case of colMeans()
#   defined by user or other packages.
# 2006-03-24
# o Created to please R CMD check.
############################################################################
