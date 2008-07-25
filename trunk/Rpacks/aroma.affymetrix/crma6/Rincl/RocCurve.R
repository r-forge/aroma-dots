setConstructorS3("RocCurve", function(roc=NULL, cuts=NULL, ...) {
  extend(BasicObject(), "RocCurve",
    roc=roc,
    cuts=cuts
  )
})

setMethodS3("as.character", "RocCurve", function(object, ...) {
  s <- sprintf("%s: ", class(object)[1]);
  s <- c(s, sprintf("AUC: %.4g", auc(object)));
  s <- c(s, sprintf("Number of cuts: %d", length(object$cuts)));
  class(s) <- "GenericSummary";
  s;
})


setMethodS3("getFpRate", "RocCurve", function(object, ...) {
  object$roc[,"fpRate"]; 
})

setMethodS3("getTpRate", "RocCurve", function(object, ...) {
  object$roc[,"tpRate"]; 
})


setMethodS3("auc", "RocCurve", function(object, ...) {
  x <- getFpRate(object);
  y <- getTpRate(object);
  ROC::trapezint(x, y, a=0, b=1);
})


setMethodS3("plot", "RocCurve", function(object, xlim=c(0,1), ylim=c(1-diff(xlim),1), lwd=2, ...) {
  plot(object$roc, type="l", xlim=xlim, ylim=ylim, lwd=lwd, ...);
})

setMethodS3("points", "RocCurve", function(object, ...) {
  points(object$roc, ...);
})

setMethodS3("lines", "RocCurve", function(object, lwd=2, ...) {
  lines(object$roc, lwd=lwd, ...);
})



############################################################################
# HISTORY:
# 2007-08-20
# o Added file caching to fitRoc2().
# 2007-08-19
# o Renamed argument 'call' to 'toCall' in fitRoc().
# 2007-04-15
# o Added scanTpAtFp().
# 2007-04-14
# o Added interpolation in findTpAtFp().
# o Removed gc() from fitRoc().
# o Added findTpAtFp() to locate TP rate for given FP rate.
# 2007-03-2x
# o Added fitRoc().
# o Created.
############################################################################
