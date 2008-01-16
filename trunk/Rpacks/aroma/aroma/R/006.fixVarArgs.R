# Add '...' to some base functions. These will later be
# turned into default functions by setMethodS3().

append <- appendVarArgs(append)
as.data.frame <- appendVarArgs(as.data.frame)
as.matrix <- appendVarArgs(as.matrix)
#log <- appendVarArgs(log)
mad <- appendVarArgs(mad)
median <- appendVarArgs(median)
power <- appendVarArgs(power)
read.table <- appendVarArgs(read.table)
var <- appendVarArgs(var)


############################################################################
# HISTORY:
# 2005-02-28
# o Created to please R CMD check.
############################################################################
