savehistory()
library(aroma.affymetrix);
library(R.graphics);
library(R.native);  # rowMedians(..., na.rm=TRUE)
source("fitRoc.R");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
shakyText <- function(x, y, labels, col=NULL, bgcol="white", adj=c(0.5,0.5), jitter=c(0.1,0.1), ...) {
  for (dx in seq(from=-jitter[1], to=+jitter[1], length=9)) {
    for (dy in seq(from=-jitter[2], to=+jitter[2], length=9)) {
      text(x, y, labels, col=bgcol, adj=adj+c(dx,dy), ...);
    }
  }
  text(x, y, labels, col=col, adj=adj, ...);
} # shakyText()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Graphical settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
imgFormat <- "screen";
imgFormat <- "eps";

mavgR <- 1;

d <- 1; by <- 0.05;

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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
case <- "optimal";

# Classification
group <- c("male"=0, "female"=1);
toCall <- 0;
pName <- names(group[group == toCall]);
nName <- names(group[group != toCall]);
tpName <- pName;
ylab <- sprintf("True-positive rate\n(correctly calling %ss %ss)", pName, pName);
xlab <- sprintf("False-positive rate\n(incorrectly calling %ss %ss)", nName, pName);




# Arrays with strong spatial artifacts (in residuals)
badArrays <- c("NA06985", "NA07055", "NA11829", "NA11830", "NA11995", "NA12004", "NA12005", "NA12006", "NA12057", "NA12234", "NA12156", "NA12236", "NA12239", "NA12264", "NA12760", "NA12762", "NA12763", "NA12892");

if (!exists("sets", mode="list")) {
  sets <- list(
   "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[1], name="CRMA"),
   "dChip,ISN,MBEI,A+B,xport" = list(col=niceCols[2], name="dChip"),
   "CNAG,FLN+GC,r60,xport" = list(col=niceCols[3], name="CNAG*"),
   "CNAT,QN,PLIER,FLN+GC,xport" = list(col=niceCols[4], name="CNAT v4")
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

# Use pool of females and males for the reference set
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
      if (refMethod == "all") {
        thetaRefs <- set$theta;
        set$thetaR <- apply(thetaRefs, MARGIN=1, FUN=median, na.rm=TRUE);
      } else if (refMethod == "adj") {
        thetaRm <- apply(set$theta[,males], MARGIN=1, FUN=median, na.rm=TRUE);
        thetaRf <- apply(set$theta[,females], MARGIN=1, FUN=median, na.rm=TRUE);
        c <- median(thetaRf, na.rm=TRUE) / median(thetaRm, na.rm=TRUE);
print(c);
        thetaRm <- thetaRm * c;
        w <- c(females=length(females), males=length(males));
        w <- w / sum(w);
        set$thetaR <- w[1] * thetaRf + w[2] * thetaRm;
      } else {
        thetaRefs <- set$theta[,references];
        set$thetaR <- apply(thetaRefs, MARGIN=1, FUN=median, na.rm=TRUE);
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
    if (sum(keep) < 1/4)
      stop(sprintf("Something is wrong. Too few loci: %.2f%%", 100*sum(keep)));
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
# Plot single-SNP ROC curves
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# One poor and one moderate SNP
snpNames <- c("SNP_A-1884486", "SNP_A-1790385");
for (snpName in snpNames) {
  jj <- indexOf(cdf, snpName);
  
  if (!Device$isOpen(fig <- fig + 1)) {
    Device$set(fig);
  
    Device$setStyle("PowerPoint");
  
    if (imgFormat == "eps") {
      imgName <- sprintf("ROC,%s,%s", case, snpName);
      filename <- sprintf("%s%s.%s", imgName, useColorsStr, imgFormat);
      eps(filename, width=6, height=6);
      par(mar=c(5.7,8.0,1.8,1.1)+0.1);
      par(cex=0.7, mgp=c(4.5,1.1,0));
    } else if (imgFormat == "screen") {
      par(mar=c(5,5,2,1)+0.1);
      par(cex=1, mgp=c(2.5,0.2,0));
    }
  
    par(cex.main=2, cex.axis=2, cex.lab=2, font.lab=2);
    par(xaxs="i", yaxs="i");
  
    names <- sapply(sets, FUN=function(set) set$name);
    
    # Create labels
    labels <- names;
  
    cols <- sapply(sets2, FUN=function(set) set$col);
    plot(NA, xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab, las=1);
    title(main=snpName);
    by <- 0.1;
    abline(h=seq(0,1,by=by), lty=3, lwd=1);
    abline(v=seq(0,1,by=by), lty=3, lwd=1);
    legend("bottomright", lty=ltys, lwd=3, col=cols, legend=labels, bg="#eeeeee", cex=1.6);
  #  legend("bottomright", pch=19, lty=ltys, lwd=4, col=cols, legend=labels, bg="#eeeeee");
  
    ltys <- rep(c(1,6), length.out=length(sets2));
    for (kk in seq(along=sets2)) {
      print(kk);
      name <- names(sets)[kk];
  
      # Get data for one SNP
      set <- sets2[[name]];
      M <- set$M;
    
      # Get data for one SNP
      Mj <- M[jj,cc];
      tj <- unname(truth[jj,cc]);
      n <- length(Mj);
  
      # Get data for one SNP
      roc <- fitRoc(tj, Mj, call=toCall, ncuts=600);
  #    plot(roc, line=TRUE, lwd=2, col=cols[kk], ltys=ltys[kk], add=TRUE);
      lines(1-roc@spec, roc@sens, lwd=3, col=cols[kk], lty=ltys[kk]);
  
  #    roc <- fitRoc(tj, Mj, call=toCall, ncuts=11);
  #    plot(roc, pch=19, cex=1.2, show.thresh=TRUE, add=TRUE);
    }
    box();
  
    if (imgFormat != "screen") {
      dev.off();
    }
  }
} # for (snpName ...)
