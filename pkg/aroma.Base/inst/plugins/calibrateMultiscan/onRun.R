onRun <- function(data, method="L1", assayNamePattern="^(.+)_(.+)$", groupBy="\\1", reports="html", imageSize=c(420,420), ..., pluginVersion="0.5", envir=parent.frame()) {
  signalFields <- c("intensity1", "intensity2");  

  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #L+# Figuring out what scans belongs to what assays
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get all 'spots' section and their assay id:s.
  spots <- getSections(data, "spots");
  totalNbrOfScans <- length(spots);

  allScanIds <- unlist(lapply(spots, FUN=getAssays));

  # Assert that duplicated data exists
  if (any(duplicated(allScanIds)))
    throw("Some 'spots' section are duplicated.");

  # Get the names of all scans
  secAssays <- getSection(data, "assays");
  allScans <- getName(secAssays, ids=allScanIds);

  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #L+# Extracting all assay names
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  allAssays <- gsub(assayNamePattern, groupBy, allScans);
  assays <- unique(allAssays);
  nbrOfAssays <- length(assays);

  #Lc# Identified ${nbrOfAssays} unique assay name(s), namely:
  #Lp# assays
  #L-#

  #L-#

  assayGroups <- list();
  for (assay in assays) {
    incl <- (assay == allAssays);
    assayGroups[[assay]] <- list(
      scanIds = allScanIds[incl],
      scans   = allScans[incl]
    )
  }

  paramFile <- "reports/.params.txt";
  mkdirs(getParent(paramFile));
  writeParam <- function(name, values=NULL, collapse="\t", append=TRUE) {
    if (is.list(values)) {
      for (kk in seq(length=length(values))) {
        writeParam(sprintf("%s/%s", name, names(values)[kk]), values[[kk]], 
                                         collapse=collapse, append=append);
      }
      return();
    }

    cat(file=paramFile, name, append=append);
    if (!is.null(values)) {
      values <- paste(values, collapse=collapse);
      cat(file=paramFile, ": ", values, sep="", append=TRUE);
    }
    cat(file=paramFile, "\n", append=TRUE);
  }

  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #L+# Write parameters needed by reporters to file '${paramFile}'
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  writeParam("# Parameters needed by reporters", append=FALSE);
  writeParam("assays", names(assayGroups));
  writeParam("scans", unlist(lapply(assayGroups, FUN=function(g) g$scans)));
  writeParam("channels", signalFields);
  #L-#

  rm(secAssays);
  rm(allScans);


  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #L+# Preparing for report creating, if requested
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!identical(reports, "none")) {
    figurePath <- "figures/";
    mkdirs(figurePath);
    assign("figurePath", figurePath, envir=envir);
    pngDevice <- findPngDevice();
    width <- 480;
    height <- 480;

    # library(RColorBrewer); x <- brewer.pal(n=12, name="Paired");
    pairedPalette <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928");

    # x <- 0:5; HclColor(h=x/max(x), l=0.5)
    scanPalette <- c("#CF3761", "#8C6E00", "#008C1D", "#008BB8", "#A341D8", "#CF3761");
    # x <- 0:12; HclColor(h=x/max(x), l=0.5)
    pairsPalette <- c("#CF3761", "#BB5100", "#9C6700", "#6E7800", "#068500", "#008E34", "#00927A", "#008DAF", "#007CD5", "#705AE1", "#B831CC", "#D2219F", "#CF3761");

    scanPalette <- pairedPalette;
    pairsPalette <- pairedPalette;
  }
  #L-#


  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #Lc# Creating a new 'assays' section
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  newSecAssays <- BaseFileAssays();
  setData(newSecAssays, data.frame(id=seq(assays), name=assays));


  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #L+# Calibrating ${nbrOfAssays} assays totalling ${totalNbrOfScans} scans
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  dprogress <- 0.5 * 95/nbrOfAssays;

  count <- 1;
  for (assay in assays) {
    assayScans <- spots[assay == allAssays];
    newAssay <- assayScans[[1]];
    nbrOfScans <- length(assayScans);
    if (nbrOfScans < 2) {
      #Lc# Skipping assay that was only scanned once: ${assay}
      cat("Skipping assay that was only scanned once: ${assay}");
      next;
    }

    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #L+# Calibrating assay '${assay}' based on ${nbrOfScans} scans
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
    for (cc in 1:length(signalFields)) {
      signalField <- signalFields[cc];
      #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      #L+# Channel ${signalField}
      #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      X <- NULL;
      for (kk in seq(length=nbrOfScans)) {
        #Lc# Extracting intensities from scan #${kk}
        scan <- assayScans[[kk]];
        Xt <- getData(scan, fields=signalField);
        Xt <- as.matrix(Xt);
        if (is.null(X))
          X <- matrix(NA, nrow=nrow(Xt), ncol=length(assayScans));
        X[,kk] <- Xt;
        rm(Xt);
      }

      #Lc# Extracted signals:
      #Ls# X

      gc();

      if ("html" %in% reports) {
        #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #L+# Creating before calibration graphs
        #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        #L+# Between scan M vs A
        imageName <- sprintf("assay-%s-Ch%d-raw-MvsA.png", assay, cc);
        pathname <- filePath(figurePath, imageName);
        pngDevice(pathname, height=imageSize[1], width=imageSize[2]);
        par(mar=c(4,4,1,1)+0.1);
        palette(pairsPalette);
        plotMvsAPairs(X);
        dev.off();
        #L-#

        #L+# Signal densities
        imageName <- sprintf("assay-%s-Ch%d-raw-density.png", assay, cc);
        pathname <- filePath(figurePath, imageName);
        pngDevice(pathname, height=0.618*imageSize[1], width=imageSize[2]);
        par(mar=c(4,4,1,1)+0.1);
        palette(scanPalette);
        plotDensities(X);
        dev.off();
        #L-#

        #L-#
      }

      #Lc# Calibrating signals (without averaging)
      X <- calibrateMultiscan(X, method=method, average=NULL);

      #L+# Fitted parameters
      fit <- attr(X, "modelFit");
      #Lc# a = ${a}
      #Lc# b = ${b}
      #Lc# adiag = ${adiag[1]}
      #Lc# converged = ${converged}
      #Lc# nbrOfIterations = ${nbrOfIterations}
      #L-#

      writeParam(sprintf("assays/%s/%s", assay, signalField), fit[c("a","b","adiag","converged","nbrOfIterations")]);
      
      if ("html" %in% reports) {
        #L+# Creating after calibration graphs

        #L+# Between scan M vs A
        imageName <- sprintf("assay-%s-Ch%d-multiscan-MvsA.png", assay, cc);
        pathname <- filePath(figurePath, imageName);
        pngDevice(pathname, height=imageSize[1], width=imageSize[2]);
        par(mar=c(4,4,1,1)+0.1);
        palette(pairsPalette);
        plotMvsAPairs(X);
        dev.off();
        #L-#

        #L+# Signal densities
        imageName <- sprintf("assay-%s-Ch%d-multiscan-density.png", assay, cc);
        pathname <- filePath(figurePath, imageName);
        pngDevice(pathname, height=0.618*imageSize[1], width=imageSize[2]);
        par(mar=c(4,4,1,1)+0.1);
        palette(scanPalette);
        plotDensities(X);
        dev.off();
        #L-#

        #L-#
      } # if ("html" %in% reports) {

      #Lc# Averaging signals (after generating graphs for report)
      X <- apply(X, MARGIN=1, FUN=median, na.rm=TRUE);
      X <- as.matrix(X);

      colnames(X) <- signalField;
      #Lc# Calibrated signals:
      #Ls# X

      #Lc# Updating intensities for first scan.
      setDataFields(newAssay, values=X);
      rm(X);

      gc();

      #L-#

      # Updating progess every channel calibrated.
      increase(progress, dprogress);
    } # for (signalField ...)
    #L-#

    #Lc# Setting parents for "merged" assay
    setParents(newSecAssays, id=count, parents=allScanIds[assay == allAssays]);

    #Lc# Setting new assay id for "merged" assay
    setAssays(newAssay, id=count);

    #Lc# Exclude all other scans.
    excl <- which(assay == allAssays)[-1];
    spots[excl] <- NA;
    rm(excl);

    count <- count + 1;
  } # for (assay ...)

  #L-#

  #Lc# Remove all "excluded" scans
  excl <- unlist(lapply(spots, FUN=identical, NA));
  spots <- spots[!excl];

  # Returning modified data
  c(list(assays=newSecAssays), spots);
}
