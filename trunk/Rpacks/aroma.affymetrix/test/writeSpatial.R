source("init.R")

verbose <- Arguments$getVerbose(TRUE);

# Find a png device that works regardless of X11 or not.
png <- System$findGraphicsDevice();
figPath <- "figures/";
mkdirs(figPath);

# Setup up the search path to ImageMagick
imageMagickConvert <- function(srcfile, destfile, format, options=NULL, ...) {
  pathname <- "C:/Program Files/ImageMagick-6.2.7-Q16/convert";
  # Escape everything
  pathname <- sprintf("\"%s\"", pathname);
  srcfile <- sprintf("\"%s\"", srcfile);
  destfile <- sprintf("\"%s\"", destfile);
  system(paste(pathname, options, srcfile, destfile));
}
options(imageConverter=imageMagickConvert);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
ds <- AffymetrixCelSet$fromFiles("chip_data/Nsp/");
print(ds);
cdf <- getCdf(ds);

# Write a spatial PNG images for CEL files
lapply(ds, function(cel) {
  print(cel);
  writeSpatial(cel, verbose=TRUE);
})

# Create thumbnails
#writeSpatial(cel, fmtstr="%s-spatial-010.png", options="-geometry 10%");
