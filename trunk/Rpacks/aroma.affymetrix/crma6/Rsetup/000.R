verbose <- Arguments$getVerbose(-10, timestamp=TRUE);
figPath <- "figures";
mkdirs(figPath);

setConstructorS3("SmartList", function(...) {
  this <- list(...);
  class(this) <- c("SmartList", class(this));
  this;
})

setMethodS3("[[<-", "SmartList", function(this, name, value) {
  oldValue <- this[[name]];
  if (identical(value, oldValue)) {
#    cat("Nothing changed!\n");
    invisible(this);
  } else {
    NextMethod("[[<-");
  }
})

setMethodS3("$<-", "SmartList", function(this, name, value) {
#  cat("$<-....\n");
  UseMethod("[[<-");
})

setMethodS3("[", "SmartList", function(this, ...) {
  class <- class(this);
  this <- NextMethod("[");
  class(this) <- class;
  this;
})

if (!exists("sets", mode="list")) {
  sets <- SmartList();
}


useColors <- FALSE;
useColors <- TRUE;
# Setup a palette
if (useColors) {
  cols <- RColorBrewer::brewer.pal(12, "Paired");
  colors <- c(CRMA=cols[6], dChip=cols[2], APT=cols[8], GTC=cols[1]);
  colors <- c(colors, setdiff(cols, colors));
  useColorsStr <- "-col";
} else {
  colors <- gray(seq(from=0, to=1, length=12));
  useColorsStr <- "";
}
palette(colors);


tpLab <- "True-positive rate\n(correctly calling males males)";
fpLab <- "False-positive rate\n(incorrectly calling males females)";

imgScales <- c("x11"=7, "eps"=6, "png"=840);
imgFormats <- names(imgScales);
force <- FALSE;
addLabels <- FALSE;
