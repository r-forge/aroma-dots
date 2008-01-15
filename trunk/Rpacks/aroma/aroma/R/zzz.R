# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
#
# Currently the following object(s) are masked from package:com.braju.sma:
# 1) put - exactly the same as in com.braju.sma.
.conflicts.OK <- TRUE


.First.lib <- function(libname, pkgname) {
  pkg <- Package(pkgname);
  assign(pkgname, pkg, pos=getPosition(pkg));
  cat(getName(pkg), " v", getVersion(pkg), " (", getDate(pkg), ")",
      " successfully loaded. See ?", pkgname, " for help.\n", sep="");
}




