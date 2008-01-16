setConstructorS3("MultiwayArray", function(x=NULL, dim=NULL, dimnames=NULL) {
  # A (nasty) workaround for the case when
  #  library(aroma) -> quit("yes") -> R -> library(R.oo2) 
  # -> Error in prod(dim) : no applicable method for "prod"
  prod <- get("prod", mode="function", pos=which(search() == "package:base"));
  
  if(is.null(x)) {
    x <- 0;
  } else {
    x <- rep(x, length.out=prod(dim));
    dim(x) <- dim;
    dimnames(x) <- dimnames;
  }
  extend(Object(x), "MultiwayArray");
})


setMethodS3("print", "MultiwayArray", function(x, ...) {
  # To please R CMD check...
  this <- x;

  tmp <- matrix(this);
  dim(tmp) <- dim(this);
  print(tmp, ...);
})

setMethodS3("apply", "MultiwayArray", function(this, MARGIN=c(), FUN, ...) {
  if (is.numeric(MARGIN))
    apply.default(this, MARGIN=MARGIN, FUN=FUN, ...)
  else
    FUN(as.vector(this), ...);
})

setMethodS3("diff", "MultiwayArray", function(x, MARGIN=c(), ...) {
  # To please R CMD check...
  this <- x;

  apply(this, MARGIN=MARGIN, FUN=diff, ...);
})

setMethodS3("mad", "MultiwayArray", function(this, MARGIN=c(), na.rm=TRUE) {
  apply(this, MARGIN=MARGIN, FUN=mad.default, na.rm=na.rm);
})

setMethodS3("mean", "MultiwayArray", function(x, MARGIN=c(), na.rm=TRUE) {
  # To please R CMD check...
  this <- x;

  apply(this, MARGIN=MARGIN, FUN=mean.default, na.rm=na.rm);
})

setMethodS3("median", "MultiwayArray", function(x, MARGIN=c(), na.rm=TRUE) {
  # To please R CMD check
  this <- x;

  apply(this, MARGIN=MARGIN, FUN=median.default, na.rm=na.rm);
})

setMethodS3("prod", "MultiwayArray", function(this, MARGIN=c(), na.rm=TRUE) {
  apply(this, MARGIN=MARGIN, FUN=prod, na.rm=na.rm);
})

setMethodS3("sum", "MultiwayArray", function(this, MARGIN=c(), na.rm=TRUE) {
  apply(this, MARGIN=MARGIN, FUN=sum, na.rm=na.rm);
})

setMethodS3("var", "MultiwayArray", function(this, MARGIN=c(), na.rm=TRUE) {
  apply(this, MARGIN=MARGIN, FUN=var.default, na.rm=na.rm);
})

setMethodS3("nbrOfDimensions", "MultiwayArray", function(this) {
  length(dim(this));
})

############################################################################
# HISTORY:
# 2002-10-24
# o Extending Object() with core value. Be careful!!!
# 2002-09-12
# o This class is "identical" to the same class in com.braju.sma.
# o Removed the usage of ObjectO2()
# 2002-01-24
# * Rewritten to make use of setMethodsS3.
# 2001-11-18
# * Added modifiers.
# 2001-11-15
# * Added sum() and prod(). Also nbrOfDimensions(),
# 2001-11-14
# * Created GSRArray.
# * Created classes MultiwayArray and Matrix.
# * Created!
############################################################################
