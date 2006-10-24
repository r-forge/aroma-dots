###########################################################################/**
# @RdocDefault writePixel
#
# @title "Creates a 1x1 bitmap image of a specific color"
#
# \description{
#  @get "title".
# }
# 
# @synopsis
#
# \arguments{
#   \item{col}{A @character string specifying the color."}
#   \item{filename}{The filename of the image.}
#   \item{path}{The path to the image.}
#   \item{ext}{The filename extension.}
#   \item{dev}{The graphical device to be used.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Return (invisibly) the pathname to the written image file.
# }
# 
# \example{\dontrun{
#   # Create a pixel image for each of the colors in the current palette
#   sapply(palette(), writePixel, path=tempdir())
# }}
#
# @author
#*/###########################################################################
setMethodS3("writePixel", "default", function(col, filename=sprintf("pixel_%s.%s", gsub("#", "", col), ext), path=NULL, ext="png", dev=png, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path);

  # Suppress warnings of 'width=1, height=1' are unlikely values in pixels.
  suppressWarnings({
    dev(pathname, width=1, height=1, bg=col);
  })
  on.exit(dev.off());
  par(mar=c(0,0,0,0));
  plot.new();

  invisible(pathname);
})

############################################################################
# HISTORY:
# 2006-10-22
# o Created.
############################################################################
