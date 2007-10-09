# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (interactive())
  savehistory();
library(aroma.affymetrix);
library(R.graphics);
library(R.native);  # rowMedians(..., na.rm=TRUE)
source("fitRoc.R");
log <- Arguments$getVerbose(-5);
timestampOn(log);


imgFormat <- "screen";
imgFormat <- "eps";
force <- FALSE;

d <- 0.06; by <- 0.01;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
getBlockAverageMap <- function(n, h=1, s=0, ...) {
  # Argument 'h':
  h <- Arguments$getDouble(h, range=c(1,1000));


  # Is h an integer?
  if (h == as.integer(h)) {
    idxs <- matrix(seq_len(n), nrow=h);
    if (n %% h != 0)
      idxs <- idxs[,-ncol(idxs)];
  } else {
    K <- ceiling(n/h);
    idxs <- matrix(TRUE, nrow=ceiling(h), ncol=K);
    steps <- (h %% 1) * ceiling(K);
    incl <- seq(from=1, to=K, length=steps);
    incl <- round(incl);
    idxs[ceiling(h), -incl] <- FALSE;
    # Shift
    if (s > 0) {
      lastRow <- idxs[ceiling(h),];
      tail <- seq(from=length(lastRow)-s+1, to=length(lastRow));
      lastRow <- c(lastRow[tail], lastRow[-tail]);
      idxs[ceiling(h),] <- lastRow;
    }
    idxs[idxs] <- seq_len(n);
    idxs[idxs == 0] <- NA;
  }

  # Skip last column in case looping
  if (n %% h != 0)
    idxs <- idxs[,-ncol(idxs)];

  # The effective 'h'
  hApprox <- sum(!is.na(idxs))/ncol(idxs);
  attr(idxs, "hApprox") <- hApprox;

  idxs;
} # getBlockAverageMap()


blockAvg <- function(X, idxs, FUN=rowMedians.matrix, ...) {
  na.rm <- (any(is.na(X)) || any(is.na(idxs)));
  dimnames <- dimnames(X);
  dimnames(X) <- NULL;
  X <- t(X);
  Z <- apply(idxs, MARGIN=2, FUN=function(jj) {
    jj <- jj[is.finite(jj)];
    Zjj <- X[,jj,drop=FALSE];
    Zjj <- FUN(Zjj, na.rm=na.rm);
    Zjj;        
  });
  Z <- t(Z);
  colnames(Z) <- dimnames[[2]];
  Z;
} # blockAvg()


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

#  cdf <- AffymetrixCdfFile$fromChipType("Mapping250K_Nsp");
#  gs <- getUnitSizes(cdf, units=units);
}

males <- which(truth[1,] == 0);
females <- which(truth[1,] == 1);


# Subset of arrays to analyze
excl <- c("NA12145");
cc <- cc0 <- which(!colnames(truth) %in% excl);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

R <- 5;

case <- "avg";


# Classification
group <- c("male"=0, "female"=1);
toCall <- 0;
pName <- names(group[group == toCall]);
nName <- names(group[group != toCall]);
tpName <- pName;
ylab <- sprintf("True-positive rate\n(correctly calling %ss %ss)", pName, pName);
xlab <- sprintf("False-positive rate\n(incorrectly calling %ss %ss)", nName, pName);

# Reference set
refSet <- "all";



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
useColors <- FALSE;
# Setup a palette
if (useColors) {
  niceCols <- RColorBrewer::brewer.pal(12, "Paired");
  niceCols <- niceCols[c(6,2)];
  useColorsStr <- "-col";
} else {
  niceCols <- c("#000000", "#666666");
  useColorsStr <- "";
}

fig <- 1;


# Arrays with strong spatial artifacts (in residuals)
badArrays <- c("NA06985", "NA07055", "NA11829", "NA11830", "NA11995", "NA12004", "NA12005", "NA12006", "NA12057", "NA12234", "NA12156", "NA12236", "NA12239", "NA12264", "NA12760", "NA12762", "NA12763", "NA12892");

if (!exists("sets", mode="list")) {
  sets <- list(
   "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[1], name="CRMA"),
   "dChip,ISN,MBEI,A+B,xport" = list(col=niceCols[2], name="dChip")
  );
}



references <- seq_len(ncol(truth));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
for (kk in seq(along=sets)) {
  name <- names(sets)[kk];
  verbose && enter(verbose, sprintf("Set #%d ('%s') of %d", 
                                                   kk, name, length(sets)));

  # Get set
  set <- sets[[name]];
  if (is.null(set$M)) {
    # Missing 'theta'?
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

    # Has 'theta'?
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

  verbose && exit(verbose);
} # for (kk in ...)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Average in non-overlapping blocks 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
keep <- (x > 2.8e6);
cc <- cc0;

fig <- 1;

# Reorder data by position
o <- order(x);
x0 <- x[o];


names <- names(sets);
names <- names[(regexpr("(CNA[GT]|h=)", names) == -1)];
for (name in names) {
  set <- sets[[name]];
  M0 <- set$M[o,];
  C0 <- truth[o,];
  n0 <- nrow(C0);

  for (h in 2:4) {
    print(c(name, h));
    newName <- sprintf("%s,h=%.3f", name, h);
    if (newName %in% names(sets))
       next;
    idxs <- getBlockAverageMap(n=n0, h=h);
    hApprox <- attr(idxs, "hApprox");
    set <- sets[[name]];
    set$h <- h;
    set$hApprox <- hApprox;
    set$x <- blockAvg(x0, idxs);
    set$M <- blockAvg(M0, idxs);
    set$truth <- blockAvg(C0, idxs);
    set$theta <- NULL;
    set$roc <- NULL;
    set$auc <- NULL;
  
    set$name <- sprintf("%s,h=%d", set$name, h);
    sets[[newName]] <- set;
  }
} # for (name in ...)


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
# Plot ROC curves
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!Device$isOpen(fig <- fig + 1)) {
  Device$set(fig);

  cols <- sapply(sets, FUN=function(set) set$col);
  ltys <- sapply(sets, FUN=function(set) set$lty);
  ltys <- rep(c(1,6), length.out=length(cols));
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
    par(mar=c(5.2,7.2,1,1.2)+0.1);
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
  labels <- c("CRMA", "dChip");
  labels <- rep(labels, times=4);
  suffix <- rep(sprintf("_%s", c("", 2:4)), each=2)
#  labels <- paste(labels, suffix, sep="");

  if (d < 0.20) {
    xs <- c(0.0545, 0.0555, 0.0168, 0.0205, 0.0162, 0.0173, 0.0036, 0.00265);
    ys <- c(0.946, 0.945, 0.979, 0.975, 0.9925, 0.9905, 0.996, 0.987);
    srt <- c(47, 49, 53, 58, 27, 35, 44, 79);
    adjy <- rep(c(-0.5,1.4), length.out=length(xs));
    for (kk in seq(along=xs)[3:4]) {
#     points(xs[kk], ys[kk], pch=19, col=cols[kk], cex=2, xpd=TRUE);
     shakyText(xs[kk], ys[kk], labels[kk], srt=srt[kk], cex=1.6, col=cols[kk], adj=c(0.5,adjy[kk]), xpd=TRUE);
    }
  } else {
    legend("bottomright", pch=19, lty=ltys, lwd=4, col=cols, legend=labels, bg="#eeeeee", cex=1.6);
  }

  if (imgFormat == "eps") {
    dev.off();
  }
}

