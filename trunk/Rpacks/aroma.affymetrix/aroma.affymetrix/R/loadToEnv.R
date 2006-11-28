###########################################################################/**
# @RdocDefault loadToEnv
#
# @title "Method to load objects to an environment"
#
# \description{
#   @get "title" for objects previously stored by @see "base::save".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "base::load".}
# }
#
# \value{
#  Returns an @environment containing all loaded objects.
# }
#
# @author
#
# \seealso{
#   See also @see "base::load".
# }
#
# \keyword{IO}
#*/###########################################################################
setMethodS3("loadToEnv", "default", function(...) {
  env <- new.env()
  base::load(..., env=env)
  env
}) # loadToEnv()


##############################################################################
# HISTORY:
# 2006-11-24
# o Created from Object.R in the R.oo package. This will probably be moved
#   to either R.oo or R.utils later.
##############################################################################
