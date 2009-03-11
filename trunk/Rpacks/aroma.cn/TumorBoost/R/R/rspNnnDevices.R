rspPngDev <- function(..., width=1, height=1, force=figForce) { 
  scale <- 640;
  width <- scale*width;
  height <- scale*height;
  pngDev(..., width=width, height=height);
}

rspEpsDev <- function(..., force=figForce) { 
  scale <- 6;
  width <- scale*width;
  height <- scale*height;
  epsDev(..., width=width, height=height, path=figPath, force=force);
}


# Default is to use PNG
figDev <- rspPngDev;
