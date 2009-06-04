onRun <- function(data, method="L1", constraint=0.05, multiassay="singleassay", referenceChannel=1, reports="html", imageSize=c(420,420), ..., pluginVersion="0.5") {
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
  #L+# Preparing for report creating, if requested
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!identical(reports, "none")) {
    figurePath <- "figures/";
    mkdirs(figurePath);
    writeParam("figurePath", figurePath);
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
  #L+# Extracting all assay names
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get all 'spots' section and their assay id:s.
  spots <- getSections(data, "spots");
  allAssayIds <- unlist(lapply(spots, FUN=getAssays));

  # Assert that duplicated data exists
  if (any(duplicated(allAssayIds)))
    throw("Some 'spots' section are duplicated.");

  # Get the names of all assays
  secAssays <- getSection(data, "assays");
  assays <- getName(secAssays, ids=allAssayIds);
  nbrOfAssays <- length(assays);

  #Lc# Identified ${nbrOfAssays} unique assay name(s), namely:
  #Lp# assays
  #L-#


  signalFields <- c("intensity1", "intensity2");  
  signalFields <- c("FCh1Median", "FCh2Median");  

  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #L+# Write parameters needed by reporters to file '${paramFile}'
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  writeParam("# Parameters needed by reporters", append=FALSE);
  writeParam("assays", assays);
  writeParam("channels", signalFields);
  #L-#



  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #L+# Normalizing ${nbrOfAssays} assays
  #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (multiassay == "singleassay") {
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #L+# Normalizing array by array
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dprogress <- 95/nbrOfAssays;

    for (spot in spots) {
      assay <- getAssays(spot)[1];
      #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      #L+# Normalizing single assay '${assay}'
      #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      #Lc# Extracting intensities
      #Lz# spot
      #Lz# signalFields
      X <- getData(spot, fields=signalFields);
      #Lc# Extracting intensities done
      X <- as.matrix(X);
 
      if ("html" %in% reports) {
        #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #L+# Creating before calibration graphs
        #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        #L+# Between scan M vs A
        imageName <- sprintf("assay-%s-raw-MvsA.png", assay);
        pathname <- filePath(figurePath, imageName);
        pngDevice(pathname, height=imageSize[1], width=imageSize[2]);
        par(mar=c(4,4,1,1)+0.1);
        palette(pairsPalette);
        plotMvsA(X);
        dev.off();
        #L-#

        #L+# Signal densities
        imageName <- sprintf("assay-%s-raw-density.png", assay);
        pathname <- filePath(figurePath, imageName);
        pngDevice(pathname, height=0.618*imageSize[1], width=imageSize[2]);
        par(mar=c(4,4,1,1)+0.1);
        palette(scanPalette);
        plotDensities(X);
        dev.off();
        #L-#

        #L-#
      }

      #Lc# Calling normalization function
      X <- normalizeAffine(X, method=method, constraint=constraint);
      colnames(X) <- signalFields;

      #L+# Fitted parameters
      fit <- attr(X, "modelFit");
      #Lc# a = ${a}
      #Lc# b = ${b}
      #Lc# converged = ${converged}
      #Lc# nbrOfIterations = ${nbrOfIterations}
      #L-#

      writeParam(sprintf("assays/%s/%s", assay, signalField), fit[c("a","b","converged","nbrOfIterations")]);

      if ("html" %in% reports) {
        #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #L+# Creating before calibration graphs
        #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        #L+# Between scan M vs A
        imageName <- sprintf("assay-%s-normalizeAffine-MvsA.png", assay);
        pathname <- filePath(figurePath, imageName);
        pngDevice(pathname, height=imageSize[1], width=imageSize[2]);
        par(mar=c(4,4,1,1)+0.1);
        palette(pairsPalette);
        plotMvsA(X);
        dev.off();
        #L-#

        #L+# Signal densities
        imageName <- sprintf("assay-%s-normalizeAffine-density.png", assay);
        pathname <- filePath(figurePath, imageName);
        pngDevice(pathname, height=0.618*imageSize[1], width=imageSize[2]);
        par(mar=c(4,4,1,1)+0.1);
        palette(scanPalette);
        plotDensities(X);
        dev.off();
        #L-#

        #L-#
      }

      #Lc# Updating intensities
      setDataFields(spot, values=X);
  
      increase(progress, dprogress);
      #L-#
    }
    #L-#
  } else if (multiassay == "multiassay") {
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #L+# Normalizing all assays together
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Allocate matrix that contain all channels
    spot <- spots[[1]];
    X <- getData(spot, fields=signalFields);
    nbrOfSpots <- nrow(X);
    X <- matrix(NA, nrow=nbrOfSpots, ncol=2*nbrOfAssays);

    # Collect data from all assays
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #L+# Collecting data from all assays
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (kk in seq(length=nbrOfAssays)) {
      spot <- spots[[kk]];
      #Lc# Assay '${getAssays(spot)[1]}'
      idx <- 2*kk - 1 + c(0,1);
      X[,idx] <- as.matrix(getData(spot, fields=signalFields));
      gc();
    }
    #L-#

    #Lc# Calling normalization function
    X <- normalizeAffine(X, method=method, constraint=constraint);

    # Collect data from all assays
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #L+# Updating intensities
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (kk in seq(length=nbrOfAssays)) {
      spot <- spots[[kk]];
      #Lc# Assay '${getAssays(spot)[1]}'
      idx <- 2*kk - 1 + c(0,1);
      setDataFields(spot, values=X[,idx], fields=signalFields);
      gc();
    }
    #L-#
    rm(X); 
    gc();

    #L-#
  } else if (multiassay == "multiassayPerChannel") {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #L+# Normalizing the channels seperately across all assays
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (nbrOfAssays < 2) {
      throw("To normalize each channel seperately across assays, at least two assays must used: ", nbrOfAssays);
    }

    # Allocate matrix that contain all channels
    spot <- spots[[1]];
    X <- getData(spot, fields=signalFields[1]);

    nbrOfSpots <- nrow(as.matrix(X));
    X <- matrix(NA, nrow=nbrOfSpots, ncol=nbrOfAssays);

    for (field in signalFields) {
      #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      #L+# Channel ${field}
      #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      #Lc# Collecting data from all assays: ${field}
      for (kk in seq(length=nbrOfAssays)) {
        spot <- spots[[kk]];
        #Lc# Assay '${getAssays(spot)[1]}'
        X[,kk] <- as.matrix(getData(spot, fields=field));
        gc();
      }

      #Lc# Calling normalization function
      X <- normalizeAffine(X, method=method, constraint=constraint);

      #Lc# Updating intensities
      for (kk in seq(length=nbrOfAssays)) {
        spot <- spots[[kk]];
        ##: Assay '${getAssays(spot)[1]}'
        setDataFields(spot, values=X[,kk,drop=FALSE], fields=field);
        gc();
      }

      #L-#
    } # for (field in ...)

    rm(X); 
    gc();
    #L-#
  } else if (multiassay == "multiassayWithReferenceChannel") {
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #L+# Normalizes the reference channel on assays together, and         \
    #L # then the other channel towards each reference channel            \
    #L # seperately
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (nbrOfAssays < 2) {
      throw("To normalize the reference channel on all assays together, at least two assays must used: ", nbrOfAssays);
    }

    dprogress <- 1/2*95/nbrOfAssays;

    # Allocate matrix that contain all channels
    spot <- spots[[1]];
    X <- getData(spot, fields=signalFields[1]);

    nbrOfSpots <- nrow(as.matrix(X));
    X <- matrix(NA, nrow=nbrOfSpots, ncol=nbrOfAssays);

    field <- signalFields[referenceChannel];
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #L+# Reference channel: ${field}
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Collect data from all assays
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #L+# Collecting data from all assays: ${field}
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (kk in seq(length=nbrOfAssays)) {
      spot <- spots[[kk]];
      #Lc# Assay '${getAssays(spot)[1]}'
      X[,kk] <- as.matrix(getData(spot, fields=field));
      gc();
    }
    #L-#

    increase(progress, 1/5*dprogress);

    #Lc# Calling normalization function
    X <- normalizeAffine(X, method=method, constraint=constraint);
    increase(progress, 3/5*dprogress);

    # Update data on all assays
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #L+# Updating intensities
    #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (kk in seq(length=nbrOfAssays)) {
      spot <- spots[[kk]];
      #Lc# Assay '${getAssays(spot)[1]}'
      setDataFields(spot, values=X[,kk,drop=FALSE], fields=field);
      gc();
    }
    #L-#

    #L-#

    rm(X); 
    gc();

    increase(progress, 1/5*dprogress);

    baselineChannel <- 1;
    if (referenceChannel == 1)
      baselineChannel <- 2;

    for (spot in spots) {
      #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      #L+# Normalizing other channel towards reference channels   \
      #L # for every single assay '${getAssays(spot)[1]}'
      #### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      #Lc# Extracting intensities
      X <- getData(spot, fields=signalFields);
      X <- as.matrix(X);
  
      #Lc# Calling normalization function
      X <- normalizeAffine(X, method=method, constraint="baseline", baselineChannel=baselineChannel);
      colnames(X) <- signalFields;

      #Lc# Updating intensities
      setDataFields(spot, values=X[,baselineChannel,drop=FALSE], fields=signalFields[baselineChannel]);

      increase(progress, dprogress);
      #L-#
    }

    #L-#
  }
  gc();

  #L-#

  # Returning modified data
  spots;
}
