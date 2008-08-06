source("init.R");

imgFormat <- c("x11", "eps", "png")[3];

fpRate <- 2 * 1/2*1/29;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit ROCs, if missing.  Plot as we go.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subset <- NULL;
J <- nrow(getData(sets[[1]]$rocData, raw=TRUE));
subset <- sample(J, size=10000);
for (unitType in knownUnitTypes[1]) {
  verbose && enter(verbose, "Estimating TP rates for each unit");
  verbose && cat(verbose, "Unit type: ", unitType);


  # Subset by unit type?
  if (unitType %in% c("snp", "cn")) {
    types <- attr(unitTypes, "map");
    rows <- whichVector(unitTypes == types[unitType]);
  } else {
    rows <- NULL;
  }

  if (!is.null(subset) && !is.null(rows)) {
    subset <- intersect(subset, 1:length(rows));
    rows <- rows[subset];
  }

  verbose && cat(verbose, "Rows:");
  verbose && str(verbose, rows);
  verbose && print(verbose, range(rows));

  imgTags <- c(unitType, "tpDensity");
  imgName <- paste(c("CRMA6", imgTags), collapse=",");

  figPathT <- filePath(figPath, "tpDensity");

  if (!devIsOpen(imgName)) {
    width <- imgScales[imgFormat] * 1;
    height <- 0.618*width;
    if (imgFormat == "x11") {
      devOpen(imgFormat, width=width, height=height, name=imgName);
	    } else {
      mkdirs(figPathT);
      imgPathname <- filePath(figPathT, sprintf("%s.%s", imgName, imgFormat));
      devOpen(imgFormat, imgPathname, width=width, height=height, name=imgName);
    }
    
    # Create empty plot
    par(mar=c(5,5,1.2,2)+0.1);
    main <- sprintf("Unit type: %s", toupper(unitType));
    ylim <- c(0,7);

    sets <- updateGraphics(sets);
    idxs <- whichVector(isRawDataSet(sets));
    isFirst <- TRUE;
    for (kk in rev(seq(along=idxs))) {
      idx <- idxs[kk];
      set <- sets[[idx]];
      verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", kk, set$name, length(idxs)));
      verbose && cat(verbose, "Name: ", set$name);
      set <- addRocData(set);
      sets[[idx]] <- set;
      rocData <- set$rocData;

      plotTpRateDensity(rocData, fpRate=fpRate, rows=rows, ylim=ylim, col=set$col, lty=set$lty, lwd=3, annotate=FALSE, add=!isFirst, verbose=verbose);

      if (isFirst) {
        stext(side=3, pos=0, cex=1, sprintf("FP rate: %.2f%%", 100*fpRate));
        n <- length(rows);
        if (is.null(rows))
          n <- dim(rocData)[1];
        stext(side=3, pos=1, cex=1, sprintf("%s units (n=%d)", toupper(unitType), n));
      }
      isFirst <- FALSE;
      verbose && exit(verbose);
    }

    devDone();
  }
  verbose && exit(verbose);
} # for (unitType ...)
gc();
