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
    # R.utils::png2()
    devices <- c(devices, png2);

    # Cairo::CairoPNG()
    if (require("Cairo")) {
      devices <- c(devices, Cairo::CairoPNG);
    }

    # grDevices::png()
    devices <- c(devices, grDevices::png);
  }

  
  System$findGraphicsDevice(devices=devices);
}, protected=TRUE)


##############################################################################
# HISTORY:
# 2007-11-25
# o We'll wait a bit making the Cairo PNG device the default *non-transparent*
#   device; the reason for this is that I still haven't figured out what the
#   all requirements are and if it is possible to use it in a "nohup" batch
#   mode without an X11 device.
# 2007-10-11
# o Created.
##############################################################################
