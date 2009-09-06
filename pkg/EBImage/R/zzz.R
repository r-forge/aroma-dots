.First.lib <- function(libname, pkgname) {
  pd <- utils::packageDescription(pkgname);
  packageStartupMessage(sprintf("Loaded %s v%s (%s). HOWEVER, %s.", pd$Package, pd$Version, pd$Date, pd$Description));
}
