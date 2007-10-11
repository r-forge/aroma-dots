setMethodS3("findPngDevice", "default", function(transparent=TRUE, ...) {
  devices <- list();
  
  if (transparent) {
    # Cairo::CairoPNG()
    if (require("Cairo")) {
      CairoPNGtrans <- function(...) {
        Cairo::CairoPNG(..., bg=NA);
      };
      devices <- c(devices, CairoPNGtrans);
    }

    # R.utils::png2()
    png2trans <- function(...) {
      png2(..., type="pngalpha");
      par(bg=NA);
      # The 'pngalpha' ghostscript device is quite slow, so to avoid
      # overloading the CPU, we add an ad hoc sleep here.
#      Sys.sleep(0.3);
    }
    devices <- c(devices, png2trans);

    # grDevices::png()
    pngtrans <- function(...) {
      png(..., bg=NA);
    }
    devices <- c(devices, pngtrans);
  } else {
    # Cairo::CairoPNG()
    if (require("Cairo")) {
      devices <- c(devices, Cairo::CairoPNG);
    }

    # R.utils::png2()
    devices <- c(devices, png2);

    # grDevices::png()
    devices <- c(devices, grDevices::png);
  }

  
  System$findGraphicsDevice(devices=devices);
}, protected=TRUE)


##############################################################################
# HISTORY:
# 2007-10-11
# o Created.
##############################################################################
