# given an input Image, transform to new RGB image based
# on list of colours given in argument 'palette' and lower and
# upper bounds given in 'lim'.
setMethodS3("rgbTransform", "Image", function(this, palette=gray.colors(256), lim=c(-Inf,Inf), outlierCol="white", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  col2rgbComposite <- function(color) {
    rgb <- col2rgb(color);
    rgb[3,]*256^2 + rgb[2,]*256 + rgb[1,];
  } # col2rgbComposite()
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'lim':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Bin the signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Binning image signals");
  # Data
  x <- this@.Data;

  # Outliers
  if (lim[1] > -Inf)
    x[x < lim[1]] <- NA;
  if (lim[2] < +Inf)
    x[x > lim[2]] <- NA;
  x[!is.finite(x)] <- NA;

  # Standardize to [0,1]
  r <- range(x, na.rm=TRUE);
  x <- (x - r[1])/(r[2]-r[1]);

  # Standardize to [0,n] where n is number of bins
  n <- length(palette);
  x <- x*n;

  # Bin to [1,n] by moving any values at the left extreme to bin one.
  x <- as.integer(x);
  x[x == 0] <- as.integer(1);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Assign colors to each bin
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Binning image data to color map");
  binValue <- col2rgbComposite(palette);
  verbose && cat(verbose, "Color map:");
  verbose && str(verbose, binValue);
  x <- binValue[x];
  x[is.na(x)] <- col2rgbComposite(outlierCol);
  verbose && str(verbose, x);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Create a new Image object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  img <- Image(x, dim(this@.Data), colormode=TrueColor);

  verbose && exit(verbose);

  img;
})


############################################################################
# HISTORY:
# 2007-01-30 /HB
# o Improved the speed even further by dirt simple binning.
# o Move to Image.rgbTransform.R.
# o Improved the speed on the internal binning in rgbTransform().
# 2007-01-12 /KS
# o TOFIX: THIS IS CURRENTLY TOO SLOW.
# o Created.  Used by getImage() in AffymetrixCelFile.
############################################################################
