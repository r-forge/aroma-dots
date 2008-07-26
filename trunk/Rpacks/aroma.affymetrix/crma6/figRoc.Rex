source("init.R");

imgFormat <- c("x11", "eps", "png")[1];

hs <- NULL;
#hs <- 2:4;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit ROCs, if missing.  Plot as we go.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ds <- c(0.10, 0.15, 0.25, 0.5, 1);
if (imgFormat == "x11") {
  ds <- ds[3];
}
for (d in ds) {
  if (0 < d && d <= 0.25) {
    by <- 0.01;
  } else if (0.25 < d && d <= 0.5) {
    by <- 0.05;
  } else {
    by <- 0.1;
  }
  for (unitType in knownUnitTypes[2]) {
    verbose && enter(verbose, "Plotting ROC");
    idxs <- seq(along=sets);
    verbose && cat(verbose, "Unit type: ", unitType);
    sets <- addSubsetDataSets(sets, subset=unitType, verbose=verbose);

    # Smooth data?
    if (!is.null(hs)) {
      sets <- addSmoothDataSets(sets, subset=unitType, hs=hs, verbose=verbose);
    }
  
    idxs <- which(isSubsetSet(sets, subset=unitType));
  #  idxs <- which(isSubsetSet(sets, subset=unitType) & isSmoothedSet(sets));
  
    dStr <- sprintf("d=%0.2f", d);
    smoothTag <- NULL;
    if (any(isSmoothedSet(sets[idxs])))
      smoothTag <- "smooth";
    imgTags <- c(unitType, smoothTag, "ROC", dStr);
    imgName <- paste(c("CRMA6", imgTags), collapse=",");
  
    if (!devIsOpen(imgName)) {
      width <- imgScales[imgFormat] * 1;
      if (imgFormat == "x11") {
        devOpen(imgFormat, width=width, height=width, name=imgName);
  	    } else {
        imgPathname <- filePath(figPath, sprintf("%s.%s", imgName, imgFormat));
        devOpen(imgFormat, imgPathname, width=width, height=width, name=imgName);
      }
      
  
      # Create empty plot
      par(mar=c(5,5,4,2)+0.1);
      par(xaxs="i", yaxs="i");
      main <- sprintf("Unit type: %s", toupper(unitType));
      plot(NA, xlim=c(0,d), ylim=c(1-d,1), xlab=fpLab, ylab=tpLab, main=main);
  
      # Add grid
      abline(h=seq(0,1,by=by), lty=3, lwd=1);
      abline(v=seq(0,1,by=by), lty=3, lwd=1);
  
      # Draw ROC curves
      sets <- updateGraphics(sets);

      for (kk in seq(along=idxs)) {
        idx <- idxs[kk];
        set <- sets[[idx]];
        verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                                kk, set$name, length(idxs)));
        verbose && cat(verbose, "Name: ", set$name);
        fit <- fit(set$rocData, force=TRUE, verbose=verbose);
        print(fit);
        lines(fit, lwd=3, col=set$col, lty=set$lty);
        verbose && exit(verbose);
      }
      
      devDone();
    }
  
    verbose && exit(verbose);
  } # for (unitType ...)
} # for (d in ds)
gc();