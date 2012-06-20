# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE;


.onAttach <- function(libname, pkgname) { 
  pkg <- Package(pkgname);
  pos <- getPosition(pkg);
  assign(pkgname, pkg, pos=pos);

  packageStartupMessage(getName(pkg), " v", getVersion(pkg), " (", 
    getDate(pkg), ") successfully loaded. See ?", pkgname, " for help.");
}

############################################################################
# HISTORY: 
# 2011-03-10
# o Created.
############################################################################
