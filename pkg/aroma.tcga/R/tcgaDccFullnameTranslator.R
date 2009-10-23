#########################################################################/**
# @set "class=character"
# @RdocMethod tcgaDccFullnameTranslator
#
# @title "Translates a TCGA DCC fullname"
#
# \description{
#   @get "title" into a comma separated format.
# }
#
# @synopsis
#
# \arguments{
#  \item{name}{A @character string.}
#  \item{...}{Not used.}
#  \item{set}{...}
# }
#
# \value{
#   Returns a @character string.
# }
#
# @author
#
# \references{
#   [1] The TCGA Network, \emph{TCGA Data Primer - Version 1.0}, July 2008.\cr
# }
#
# \seealso{
#   @seemethod "parseTcgaDccArchiveName".
# }
#
# @keyword "file"
#*/######################################################################### 
setMethodS3("tcgaDccFullnameTranslator", "character", function(name, ..., set=NULL) {
  path <- getPath(set);
  name <- basename(path);
  parts <- parseTcgaDccArchiveName(name);
  fields <- c("domain", "tumorType", "platform", "version");
  paste(parts[fields], collapse=",");
}) # tcgaDccFullnameTranslator()


############################################################################
# HISTORY:
# 2009-06-07
# o Created.
############################################################################
