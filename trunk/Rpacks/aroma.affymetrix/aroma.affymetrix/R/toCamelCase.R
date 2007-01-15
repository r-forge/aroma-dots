###########################################################################/**
# @RdocDefault toCamelCase
#
# @title "Converts a string of words into a merged camel-cased word"
#
# \description{
#   @get "title", e.g. "a short black" is converted to "aShortBlack".
# }
#
# @synopsis
#
# \arguments{
#  \item{s}{A @character string.}
#  \item{capitalize}{If @TRUE, the first letter will be in upper case,
#    otherwise it will be in lower case.}
#  \item{split}{A pattern used to identify words.  See @see "base::strsplit"
#    for more details.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @see "R.utils::capitalize" and @see "R.utils::decapitalize".
#   @see "base::chartr".
# }
#
# @keyword programming
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("toCamelCase", "default", function(s, capitalize=FALSE, split="[ \t]+", ...) {
  s <- strsplit(s, split=split);
  s <- lapply(s, FUN=function(s) {
    paste(capitalize(tolower(s)), collapse="")
  });
  s <- unlist(s);
  if (!capitalize)
    s <- decapitalize(s);
  s;
}, private=TRUE)

############################################################################
# HISTORY:
# 2007-01-14
# o Added Rdoc comments.
# 2006-09-15
# o Created. Will probably end up in R.utils some day.
############################################################################
