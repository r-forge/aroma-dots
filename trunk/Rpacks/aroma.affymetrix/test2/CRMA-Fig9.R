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

d <- 1; by <- 0.05
d <- 0.11; by <- 0.01;
d <- 0.09; by <- 0.01;

# Setup a palette
useColors <- FALSE;
# Setup a palette
if (useColors) {
  niceCols <- RColorBrewer::brewer.pal(12, "Paired");
  niceCols <- niceCols[c(6,5,10)];
  useColorsStr <- "-col";
} else {
  niceCols <- c("#000000", "#666666", "#999999");
  useColorsStr <- "";
}

fig <- 1;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
case <- "reference";

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
   "ACC,+300,RMA,A+B,FLN,-X;all" = list(col=niceCols[1], name="all"),
   "ACC,+300,RMA,A+B,FLN,-X;adj" = list(col=niceCols[2], name="adjusted"),
   "ACC,+300,RMA,A+B,FLN,-X;diploid" = list(col=niceCols[3], name="diploid")
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
    if (nok < 1/4)
      stop(sprintf("Something is wrong. Call rate is only %.2f%%", 100*nok));
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot ROC curves
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!Device$isOpen(fig <- fig + 1)) {
  Device$set(fig);

  cols <- sapply(sets, FUN=function(set) set$col);
  ltys <- sapply(sets, FUN=function(set) set$lty);
  names <- sapply(sets, FUN=function(set) set$name);
  
  # Create labels
  labels <- names;
  #labels <- splitByCommonTails(labels)[,2];
  for (kk in seq(along=sets)) {
    labels[kk] <- sprintf("[%.2f%%] %s", 100*aucs[kk], labels[kk]);
    if (imgFormat == "screen")
      labels[kk] <- sprintf("%s (%.2f%%)", labels[kk], 100*callRates[kk]);
  }
  print(labels);
  
  xlim <- c(0,d); ylim=c(1-d,1);
  dStr <- sprintf("%0.3f", d);

  Device$setStyle("PowerPoint");
 
  if (imgFormat == "eps") {
    dStr2 <- gsub("[.]", "_", dStr);
    imgName <- sprintf("ROC,%s,h%d,d%s", case, mavgR, dStr2);
    filename <- sprintf("%s%s.eps", useColorsStr, imgName);
    eps(filename, width=5, height=5);
    par(mar=c(5.2,7.2,1,1)+0.1);
    par(cex=0.7, mgp=c(4.1,0.9,0));
  } else if (imgFormat == "screen") {
    par(mar=c(5,5,2,1)+0.1);
    par(cex=1, mgp=c(2.5,0.2,0));
  }

  # Create empty plot
  par(cex.axis=1.7, cex.lab=1.7, font.lab=2);
  par(xaxs="i", yaxs="i");
  plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, las=1);
 
  # Add grid
  abline(h=seq(0,1,by=by), lty=3, lwd=1);
  abline(v=seq(0,1,by=by), lty=3, lwd=1);

  # Draw ROC curves
  for (kk in seq(along=sets)) {
    print(labels[kk]);
    set <- sets[[kk]];
    roc <- set$roc;
    plot(roc, line=TRUE, lwd=4, lty=ltys[kk], col=cols[kk], add=TRUE);
  } # for (kk in ...)

  # Add labels/legend
  labels <- sapply(sets, FUN=function(set) set$name);
  if (d < 0.20) {
    xs <- c(0.0535, 0.0555, 0.0645);
    ys <- c(0.945, 0.945, 0.933);
    srt <- rep(45, 3);
    adjy <- c(-0.6,1.4,-0.5);
    for (kk in seq(along=xs)) {
#      points(xs[kk], ys[kk], pch=19, col=cols[kk]);
     shakyText(xs[kk], ys[kk], labels[kk], srt=srt[kk], cex=1.8, col=cols[kk], adj=c(0.5,adjy[kk]));
    }
  } else {
    legend("bottomright", pch=19, lty=ltys, lwd=4, col=cols, legend=labels, bg="#eeeeee", cex=1.6);
  }

  if (imgFormat == "eps") {
    dev.off();
  }
}

