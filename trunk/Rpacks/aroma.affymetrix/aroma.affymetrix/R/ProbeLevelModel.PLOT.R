setMethodS3("plotMvsPosition", "ProbeLevelModel", function(this, sample, ..., annotate=TRUE) {
  ces <- getChipEffects(this);
  ce <- getFile(ces, sample);
  cesAvg <- getAverageFile(ces);
  res <- plotMvsPosition(ce, reference=cesAvg, ..., annotate=annotate);  
  if (annotate) {
    stext(getLabel(this), side=1, pos=1, line=-1, cex=0.7, col="darkgrey");
  }

  this$lastPlotData <- res;
  this$lastPlotSample <- sample;

  invisible(res);
}, private=TRUE)


setMethodS3("highlight", "ProbeLevelModel", function(this, ...) {
  ces <- getChipEffects(this);
  ce <- getFile(ces, this$lastPlotSample);
  highlight(ce, ...);
}, private=TRUE)


############################################################################
# HISTORY:
# 2006-08-28
# o Added plotMvsPosition().
############################################################################
