#########################################################################/**
# @set "class=character"
# @RdocMethod parseTcgaDccArchiveName
#
# @title "Parses a TCGA DCC archive name"
#
# \description{
#   @get "title" into its components.
# }
#
# @synopsis
#
# \arguments{
#  \item{name}{A @character string.}
#  \item{unlist}{A @logical.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @character @vector of length N=6.
# }
#
# @author
#
# \examples{
#   name <- "broad.mit.edu_GBM.HG-U133_Plus_2.1.0.0"
#   parts <- parseTcgaDccArchiveName(name)
# }
#
# \references{
#   [1] The TCGA Network, \emph{TCGA Data Primer - Version 1.0}, July 2008.\cr
# }
#
# @keyword "file"
#*/######################################################################### 
setMethodS3("parseTcgaDccArchiveName", "character", function(name, unlist=FALSE, ...) {
  pattern <- "^([a-z.]+)_([A-Z]+)[.]([^.]+)[.]([0-9]+)[.]([0-9]+)[.]([0-9]+).*$";
#  pattern <- toAsciiRegExprPattern(pattern);

  parts <- list();
  parts$domain <- gsub(pattern, "\\1", name);
  parts$tumorType <- gsub(pattern, "\\2", name);
  parts$platform <- gsub(pattern, "\\3", name);
  parts$archiveSerialIndex <- gsub(pattern, "\\4", name);
  parts$revision <- gsub(pattern, "\\5", name);
  parts$series <- gsub(pattern, "\\6", name);

  fields <- c("archiveSerialIndex", "revision", "series");
  version <- paste(parts[fields], collapse=".");
  parts$version <- version;

  if (unlist) {
    parts <- unlist(parts);
  }

  parts;
})


############################################################################
# HISTORY:
# 2009-06-07
# o Created.
############################################################################
