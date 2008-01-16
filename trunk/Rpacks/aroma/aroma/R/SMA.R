#########################################################################/**
# @RdocClass SMA
#
# @title "The SMA class"
#
# \description{
#  @classhierarchy
#
#   The SMA class is a static class that provided methods to convert
#   data object in aroma into sma structures. This could be
#   useful if you would like to use methods in sma that are not (yet)
#   implemented in the aroma package.
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/#########################################################################
setConstructorS3("SMA", function() {
  extend(Object(), "SMA")
}, static=TRUE);




#########################################################################/**
# @RdocMethod as.RG
#
# @title "Converts a aroma object into an object of the RG structure"
#
# @synopsis
#
# \description{
#   @get "title", which is used by most of the function in the sma package.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("as.RG", "SMA", function(this, obj) {
  if (is.null(obj))
    throw("Argument 'obj' is NULL.");

  if (inherits(obj, "MAData")) {
    obj <- as.RGData(obj);
  }
  
  if (inherits(obj, "RawData")) {
    list(R=obj$R, G=obj$G, Rb=obj$Rb, Gb=obj$Gb);
  } else if (inherits(obj, "RGData")) {
    zeros <- matrix(0, nrow=nrow(obj$R), ncol=ncol(obj$R));
    list(R=obj$R, G=obj$G, Rb=zeros, Gb=zeros);
  } else {
    throw("Can not convert to the RG data structure. Unknown data type: ", data.class(obj));
  }
}, static=TRUE);



#########################################################################/**
# @RdocMethod as.layout
#
# @title "Converts a aroma object into an object of the layout structure"
#
# @synopsis
#
# \description{
#   @get "title", which is used by most of the functions in the sma package.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("as.layout", "SMA", function(this, obj) {
  if (is.null(obj))
    throw("Argument 'obj' is NULL.");

  if (inherits(obj, "MicroarrayData")) {
    obj <- getLayout(obj);
    if (is.null(obj))
      throw("The given object does not specify an Layout object.");
  }

  if (!inherits(obj, "Layout")) 
    throw("Can not convert to layout. Argument is of unknown type: ", data.class(obj));
  
  list(nspot.r=obj$nspot.r, nspot.c=obj$nspot.c, ngrid.r=obj$ngrid.r, ngrid.c=obj$ngrid.c);
}, static=TRUE);


setMethodS3("loadData", "SMA", function(this, names=NULL) {
  totalSize <- 0;
  require(sma);
  data(MouseArray); # Approximately 10Mb
  if (!is.null(names)) {
    # Keep only the variables specified by 'names'
    loadedNames <- c("mouse.data", "mouse.gnames", "mouse.lratio", "mouse.setup", "mouse.t2", "mouse1", "mouse2", "mouse3", "mouse4", "mouse5", "mouse6");
    residues <- setdiff(loadedNames, names);
    SMA$unloadData(residues);
  }
  # Calculate the total memory loaded.
  for (name in names) {
    if (exists(name))
      totalSize <- totalSize + object.size(get(name));
  }

  # Check if gco() should be ran. The rational for doing this here is that
  # loadData() is mostly used in Rd examples and when doing R CMD check on
  # the package the gco() is never called and R CMD check will quickly run
  # into memory problems. /HB 2002-06-24
  OST <- .Platform$OS.type;
  if (OST == "windows") {
    memory.free <- memory.limit() - memory.size();
  } else {
    memory.free <- Inf;
  }
  if (memory.free < 5e6) {
    warning("Running low of memory. Calling garbage collector gc().");
    gc();
  }
  
  invisible(totalSize);
}, static=TRUE)


setMethodS3("unloadData", "SMA", function(this, names=NULL, envir=.GlobalEnv) {
  # Saves about 10Mb if unloading everything.
  totalSize <- 0;
  if (is.null(names))
    names <- c("mouse.data", "mouse.gnames", "mouse.lratio", "mouse.setup", "mouse.t2", "mouse1", "mouse2", "mouse3", "mouse4", "mouse5", "mouse6");
  for (name in names) {
    if (exists(name, envir=envir)) {
      totalSize <- totalSize + object.size(get(name, envir=envir));
      rm(list=name, envir=envir);
    }
  }
  gc();
  invisible(totalSize);
}, static=TRUE)



############################################################################
# HISTORY:
# 2002-09-13
# o BUG FIX: memory.limit() and memory.size() do only exist on Windows so
#   for now the automatic calling of gco() only works on Windows.
# 2002-06-24
# * SMA$loadData() now calls gco() if [R] is running low of memory. The
#   rational for doing this here is that loadData() is mostly used in Rd
#   examples and when doing R CMD check on the package the gco() is never
#   called and R CMD check will quickly run into memory problems.
# 2002-05-06
# * Added SMA$loadData() and SMA$unloadata(). The SMA data set takes about
#   10Mb of memory, so it is is wise to unload it when not used. Note that
#   it is not possible to call it SMA$data() because then the generic
#   function will make data(...) stop working.
# 2002-02-27
# * Rewritten to make use of setMethodS3's.
# 2001-08-08
# * Created!
############################################################################

