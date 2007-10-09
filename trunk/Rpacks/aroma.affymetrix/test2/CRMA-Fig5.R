savehistory()
library(aroma.affymetrix);
library(R.graphics);
library(R.native);  # rowMedians(..., na.rm=TRUE)
source("fitRoc.R");
source("plotCumulativeDensity.list.R");

imgFormat <- "screen";
imgFormat <- "eps";

# Smoothing parameter (zero == no smoothing)
sd <- 1e3;  # in units of 1bp
sd <- 0;

mavgR <- 1;

  d <- 1; by <- 0.05
  d <- 0.50; by <- 0.05;
  d <- 0.25; by <- 0.05;
  d <- 0.15; by <- 0.01;
  d <- 0.12; by <- 0.01;
  d <- 0.10; by <- 0.01;

case <- "optimal";

# Classification
group <- c("male"=0, "female"=1);
toCall <- 0;
pName <- names(group[group == toCall]);
nName <- names(group[group != toCall]);
tpName <- pName;
ylab <- sprintf("True-positive rate\n(correctly calling %ss %ss)", pName, pName);
xlab <- sprintf("False-positive rate\n(incorrectly calling %ss %ss)", nName, pName);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

useColors <- FALSE;
# Setup a palette
if (useColors) {
  niceCols <- RColorBrewer::brewer.pal(12, "Paired");
  niceCols <- niceCols[c(6,2,8,4)];
  useColorsStr <- "-col";
} else {
  niceCols <- c("#000000", "#444444", "#888888", "#aaaaaa");
  useColorsStr <- "";
}

fig <- 1;


# Arrays with strong spatial artifacts (in residuals)
badArrays <- c("NA06985", "NA07055", "NA11829", "NA11830", "NA11995", "NA12004", "NA12005", "NA12006", "NA12057", "NA12234", "NA12156", "NA12236", "NA12239", "NA12264", "NA12760", "NA12762", "NA12763", "NA12892");

if (!exists("sets", mode="list")) {
  sets <- list(
   "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[1], name="CRMA"),
   "dChip,ISN,MBEI,A+B,xport" = list(col=niceCols[2], name="dChip"),
   "CNAG,FLN+GC,r60,xport" = list(col=niceCols[3], name="CNAG*"),
   "CNAT,QN,PLIER,FLN+GC,xport" = list(col=niceCols[4], name="CNAT")
  );
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dataSetName <- "Affymetrix_2006-HapMap270.CEU.founders";
dataPath <- filePath("data", dataSetName);

if (!exists("truth")) {
  filename <- "etc.RData";
  pathname <- Arguments$getReadablePathname(filename, path=dataPath, 
                                                     mustExist=TRUE);
  attachLocally(loadObject(pathname));

  filename <- "truth.RData";
  pathname <- Arguments$getReadablePathname(filename, path=dataPath, 
                                                     mustExist=TRUE);
  attachLocally(loadObject(pathname));

  cdf <- AffymetrixCdfFile$fromChipType("Mapping250K_Nsp");
#  gs <- getUnitSizes(cdf, units=units);
}

males <- which(truth[1,] == 0);
females <- which(truth[1,] == 1);

# Use pool of females and males as a reference
references <- seq_len(ncol(truth));


for (kk in seq(along=sets)) {
  name <- names(sets)[kk];
  set <- sets[[name]];
  if (is.null(set$M)) {
    if (is.null(set$theta)) {
      # Load thetas for chromosome X
      dataName <- gsub(";.*$", "", name);
      filename <- sprintf("%s.RData", dataName);
      pathname <- Arguments$getReadablePathname(filename, path=dataPath, 
                                                             mustExist=TRUE);
      print(filename);
      data <- loadObject(pathname);
      if ("theta" %in% names(data))
        set$theta <- data$theta;
      if ("M" %in% names(data))
        set$M <- data$M;
      rm(data);
    }

    if (!is.null(set$theta)) {
      # Calculate average theta, and (A,M).
      refMethod <- gsub("^.*;", "", name);
      print(refMethod);
      if (refMethod == "all") {
        print("ALL");
        thetaRefs <- set$theta;
        set$thetaR <- rowMedians.matrix(thetaRefs, na.rm=TRUE);
      } else if (refMethod == "diploid") {
        print("ALL");
        thetaRefs <- set$theta[,females];
        set$thetaR <- rowMedians.matrix(thetaRefs, na.rm=TRUE);
      } else if (refMethod == "adj") {
        print("ADJUSTED");
        thetaRm <- rowMedians.matrix(set$theta[,males], na.rm=TRUE);
        thetaRf <- rowMedians.matrix(set$theta[,females], na.rm=TRUE);
        c <- median(thetaRf, na.rm=TRUE) / median(thetaRm, na.rm=TRUE);
print(c);
        thetaRm <- thetaRm * c;
        w <- c(females=length(females), males=length(males));
        w <- w / sum(w);
        set$thetaR <- w[1] * thetaRf + w[2] * thetaRm;
      } else if (refMethod == "mean") {
        print("MEAN; ALL or DIPLOID");
        thetaRefs <- set$theta[,references];
        set$thetaR <- rowMeans(thetaRefs, na.rm=TRUE);
      } else {
        print("ALL or DIPLOID");
        thetaRefs <- set$theta[,references];
        set$thetaR <- rowMedians.matrix(thetaRefs, na.rm=TRUE);
      }
      str(thetaRefs);

      set$M <- log2(set$theta/set$thetaR); 
      set$A <- log2(set$theta*set$thetaR)/2;
    } else if (!is.null(set$M)) {
      if (name == "CNAT,QN,PLIER,FLN+GC,xport") {
        # Make sure order of samples are the same as in 'truth'
        df <- read.table("CNAT/foo.xls", header=TRUE, sep="\t");
        expId <- sprintf("%s.cn", df[,1]);
        sampleNames <- gsub("[.]cel$", "", df[,2]);
        cc <- match(colnames(set$M), expId);
        stopifnot(all(is.finite(cc)));
        sampleNames <- sampleNames[cc];
        cc <- match(sampleNames, colnames(truth));
        stopifnot(all(is.finite(cc)));
        set$M <- set$M[,cc];
        rm(df, expId, sampleNames, cc);
      }
    }

    # Assert that we at least got an 'M' field
    stopifnot("M" %in% names(set));

    sets[[name]] <- set;
  }
}




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup graphical annotations
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
for (name in names(sets)) {
  set <- sets[[name]];

  if (is.null(set$col))
    set$col <- kk;
  if (is.null(set$lty))
    set$lty <- 1;
  if (is.null(set$name))
    set$name <- name;

  sets[[name]] <- set;
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit ROC curves
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Subset of arrays to analyze
cc <- seq_len(ncol(set$M));
#cc <- which(!colnames(set$M) %in% badArrays)
# Exclude one sample that behaves poorly for dChip
#cc <- setdiff(cc, males[21]);
#cc <- setdiff(cc, females[1]);
#cc <- setdiff(cc, 11);
cc <- setdiff(cc, 32);

for (name in names(sets)) {
  print(name);

  set <- sets[[name]];

  if (is.null(set$roc)) {
    # Get data
    M <- set$M;
    A <- set$A;

    xx <- set$x;
    if (is.null(xx))
      xx <- x;

    C <- set$truth;
    if (is.null(C))
      C <- truth;

    # Exclude pseudo copy-neutral regions
    keep <- (xx > 2.8e6);
    xx <- xx[keep];
    M <- M[keep,cc,drop=FALSE];
    C <- C[keep,cc,drop=FALSE];

    # Vectorize
    M <- as.vector(M);
    C <- as.vector(C);
    
    # Keep only finite data points
    ok <- (is.finite(M) & is.finite(C));
    n <- length(M);
    nok <- sum(ok);
    callRate <- nok/n;
    cat(sprintf("Keeping %d out of %d (%.2f%%) data points \n", 
                                               nok, n, 100*callRate));

    M <- M[ok];
    C <- C[ok];
  
    # fit ROC curve
    fit <- fitRoc(C, M, call=toCall);

    set$nbrOfSnps <- n;
    set$callRate <- callRate;
    set$roc <- fit;
    rm(fit, ok, M, A);
  }

  if (is.null(set$auc))
    set$auc <- AUC(set$roc);

  sets[[name]] <- set;
}
aucs <- sapply(sets, FUN=function(set) set$auc);

# Reorder methods by AUC
o <- order(aucs, decreasing=TRUE);

aucs <- aucs[o];
sets <- sets[o];

callRates <- sapply(sets, FUN=function(set) set$callRate);



# Normalize to females
if (!exists("sets2", mode="list") || !identical(names(sets2), names(sets))) {
  sets2 <- lapply(sets, FUN=function(set) {
    M <- set$M;
    Mr <- apply(M[,females], MARGIN=1, FUN=median, na.rm=TRUE);
    M <- M - Mr;
    set$M <- M;
    set;
  })
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Density of TP rates across SNPs at fixed FP rate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!Device$isOpen(fig <- fig + 1)) {
  Device$set(fig, height=4.5);
  C <- truth[,cc];
  nbrOfFemales <- sum(truth[1,cc]);
  fpStep <- 1/nbrOfFemales;

  cols <- sapply(sets, FUN=function(set) set$col);

  fpRate <- fpStep/2;
 
  print(fpRate);

  filename <- sprintf("TPatFP%.5f.RData", fpRate);
  pathname <- filename;
  print(pathname);
  if (isFile(pathname)) {
    D <- loadObject(pathname);
  } else {
    D <- matrix(NA, nrow=nrow(C), ncol=length(sets)); 
    colnames(D) <- sapply(sets, FUN=function(set) set$name);
    attr(D, "fpRate") <- fpRate;
    
    for (mm in seq_len(ncol(D))) {
      print(mm);
      set <- sets[[mm]];
      M <- set$M[,cc];
      D[,mm] <- fitTpDensity(C, M, fpRate=fpRate, call=toCall, d=D[,mm]); 
    }
  
    saveObject(D, pathname);
  }
  
  if (imgFormat == "eps") {
    imgName <- sprintf("ROC,%s,TPdistribution,fpRate%.4f", case, fpRate);
    imgName <- gsub("[.]", "_", imgName);
    filename <- sprintf("%s%s.%s", imgName, useColorsStr, imgFormat);
    eps(filename, width=6, height=0.618*6);
    par(mar=c(4.0,4.5,1,1)+0.1);
    par(cex.axis=1.7, cex.lab=1.7, font.lab=2);
    par(cex=0.7, mgp=c(2.7,0.9,0));
  } else if (imgFormat == "screen") {
    par(mar=c(5,5,2,1)+0.1);
    par(cex=1, mgp=c(2.5,0.2,0));
  }
  
  
  xlabTp <- "TP rates (%)";
  plotDensity(100*D, from=0, to=100, xlab=xlabTp, lwd=2, col=cols, lty=c(1,6));

  # Record FP rate
  stext(side=4, pos=0, sprintf("FP rates %.2f%%", 100*fpRate), cex=0.7);

  # Add labels/legend
  labels <- sapply(sets, FUN=function(set) set$name);
  if (d < 0.20) {
    xs <- c(95, 93.5, 93.2, 91);
    ys <- c(0.05, 0.0435, 0.0373, 0.030)-0.002;
    adjx <- rep(1, 4);
    for (kk in 1:4) {
#      points(xs[kk], ys[kk], pch=19, col=cols[kk]);
      text(xs[kk], ys[kk], labels[kk], cex=1.2, col=cols[kk], adj=c(adjx[kk], 0.5));
    }
  } else {
    legend("bottomright", pch=19, lty=ltys, lwd=4, col=cols, legend=labels, bg="#eeeeee", cex=1.8);
  }
    
  if (imgFormat != "screen") {
    dev.off();
  }
}
