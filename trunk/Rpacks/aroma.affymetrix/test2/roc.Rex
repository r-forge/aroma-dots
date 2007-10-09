savehistory()
library(aroma.affymetrix);
library(R.graphics);
library(R.native);  # rowMedians(..., na.rm=TRUE)
source("fitRoc.R");
source("plotCumulativeDensity.list.R");

imgFormat <- "screen";
#imgFormat <- "eps";

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
#  d <- 0.06; by <- 0.01;
#  d <- 0.05; by <- 0.01;
#  d <- 0.025; by <- 0.001;
#  d <- 0.01; by <- 0.001;

case <- "optimal";
#case <- "steps";
#case <- "allVdiploid";

# Classification
group <- c("male"=0, "female"=1);
toCall <- 0;
pName <- names(group[group == toCall]);
nName <- names(group[group != toCall]);
tpName <- pName;
ylab <- sprintf("True-positive rate\n(correctly calling %ss %ss)", pName, pName);
xlab <- sprintf("False-positive rate\n(incorrectly calling %ss %ss)", nName, pName);

# Reference set
refSets <- c("all", "females", "males");
refSet <- "all";
#refSet <- "females";



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup a palette
niceCols <- RColorBrewer::brewer.pal(12, "Paired");

fig <- 1;


# Arrays with strong spatial artifacts (in residuals)
badArrays <- c("NA06985", "NA07055", "NA11829", "NA11830", "NA11995", "NA12004", "NA12005", "NA12006", "NA12057", "NA12234", "NA12156", "NA12236", "NA12239", "NA12264", "NA12760", "NA12762", "NA12763", "NA12892");

if (!exists("sets", mode="list")) {
  if (case == "one") {
    sets <- list(
     "QN,RMA,A+B" = list(col=niceCols[6], name="QN,RMA,A+B"),
    );
  } else if (case == "optimal") {
    sets <- list(
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[6], name="CRMA"),
#     "ACC,+300,RMA,A+B,FLN,-X,SN" = list(col=niceCols[5], name="CRMA,SN"),
#     "ACC,+300,RMA,A+B,FLN,-X,GCN,-X" = list(col=niceCols[5], name="CRMA,GC"),
#     "dChip,ISN,MBEI,A+B,00,xport" = list(col=niceCols[1], name="dChip,no-outliers"),
     "dChip,ISN,MBEI,A+B,xport" = list(col=niceCols[2], name="dChip"),
#     "dChip,QN,MBEI,A+B,xport" = list(col=niceCols[1], name="dChip,QN"),
     "CNAG,FLN+GC,r60,xport" = list(col=niceCols[8], name="CNAG*"),
     "CNAT,QN,PLIER,FLN+GC,xport" = list(col=niceCols[4], name="CNAT v4")
    );
  } else if (case == "MTV") {
    sets <- list(
#     "ACC,+350,RMA,A+B,FLN,-X" = list(col=niceCols[3], name="CRMA,+350"),
#     "ACC,+250,RMA,A+B,FLN,-X" = list(col=niceCols[4], name="CRMA,+250"),
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[5], name="CRMA"),
     "ACC,+300,RMA,A+B,FLN,-X,T" = list(col=niceCols[6], name="CRMA,T"),
     "ACC,+300,RMA,A+B,FLN,-X,V" = list(col=niceCols[7], name="CRMA,V")
    );
  } else if (case == "DCQN") {
    sets <- list(
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[5], name="CRMA"),
     "ACC,QN,pm,+300,RMA,A+B,FLN,-X" = list(col=niceCols[3], name="CRMA,QN"),
     "ACC,DCQN,pm,+300,RMA,A+B,FLN,-X" = list(col=niceCols[4], name="CRMA,DCQN"),
     "ACC,DCQN,pm,-X,+300,RMA,A+B,FLN,-X" = list(col=niceCols[7], name="CRMA,DCQN,-X"),
     "ACC,DCQN,pm,-X,robust,+300,RMA,A+B,FLN,-X" = list(col=niceCols[8], name="CRMA,DCQN,-X,robust")
    );

  } else if (case == "optimalMT") {
    sets <- list(
     "CNAG,FLN+GC,r60,xport" = list(col=niceCols[8], name="CNAG*"),
     "CNAT,QN,PLIER,FLN+GC,xport" = list(col=niceCols[4], name="CNAT v4"),
     "dChip,ISN,MBEI,A+B,xport" = list(col=niceCols[2], name="dChip"),
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[5], name="CRMA"),
     "ACC,+300,RMA,A+B,FLN,-X,T" = list(col=niceCols[6], name="CRMA,T")
    );
  } else if (case == "optimalT") {
    sets <- list(
     "CNAG,FLN+GC,r60,xport" = list(col=niceCols[8], name="CNAG*"),
     "CNAT,QN,PLIER,FLN+GC,xport" = list(col=niceCols[4], name="CNAT v4"),
     "dChip,ISN,MBEI,A+B,xport" = list(col=niceCols[2], name="dChip"),
     "ACC,+300,RMA,A+B,FLN,-X,T" = list(col=niceCols[6], name="CRMA,T")
    );
  } else if (case == "optimalT") {
    sets <- list(
     "CNAG,FLN+GC,r60,xport" = list(col=niceCols[8], name="CNAG*"),
     "CNAT,QN,PLIER,FLN+GC,xport" = list(col=niceCols[4], name="CNAT v4"),
     "dChip,ISN,MBEI,A+B,xport" = list(col=niceCols[2], name="dChip"),
     "ACC,+300,RMA,A+B,FLN,-X,T" = list(col=niceCols[6], name="CRMA,T"),
     "ACC,+300,RMA,A+B,FLN,-X,V" = list(col=niceCols[7], name="CRMA,V")
    );
  } else if (case == "MvT") {
    sets <- list(
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[2], name="CRMA"),
     "ACC,+300,RMA,A+B,FLN,-X,T" = list(col=niceCols[6], name="CRMA,T")
    );
  } else if (case == "FLN") {
    sets <- list(
     "ACC,+300,RMA,A+B" = list(col=niceCols[2], name="ACC,RMA"),
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[6], name="ACC,RMA,FLN"),
    );
  } else if (case == "steps") {
    sets <- list(
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[6], name="ACC+FLN"),
     "ACC,+300,RMA,A+B" = list(col=niceCols[4], name="ACC"),
     "QN,RMA,A+B" = list(col=niceCols[10], name="QN")
    );
  } else if (case == "allVdiploid") {
    sets <- list(
     "ACC,+300,RMA,A+B,FLN,-X;all" = list(col=niceCols[6], name="CRMA,all"),
     "ACC,+300,RMA,A+B,FLN,-X;mean" = list(col=niceCols[7], name="CRMA,mean"),
     "ACC,+300,RMA,A+B,FLN,-X;adj" = list(col=niceCols[5], name="CRMA,adj"),
     "ACC,+300,RMA,A+B,FLN,-X;diploid" = list(col=niceCols[10], name="CRMA,diploid")
#     "dChip,ISN,MBEI,A+B,xport;adj" = list(col=niceCols[1], name="dChip,adj"),
#     "dChip,ISN,MBEI,A+B,xport;diploid" = list(col=niceCols[2], name="dChip,diploid"),
    );
  } else if (case == "CQN") {
    sets <- list(
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[6], name="ACC,RMA,FLN"),
     "ACC,+300,RMA,A+B,FLN,-X,DCQN,-X,v0" = list(col=niceCols[5], name="ACC,RMA,FLN,QN,v0"),
     "ACC,+300,RMA,A+B,FLN,-X,DCQN,-X" = list(col=niceCols[4], name="ACC,RMA,FLN,QN")
    );
  } else if (case == "+-") {
    sets <- list(
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[6], name="ACC,RMA,FLN"),
     "ACC,+300,RMA,+-,A+B,FLN,-X,CEGM" = list(col=niceCols[5], name="ACC,RMA,+-,FLN,CEGM")
    );
  } else if (case == "PM-MM") {
    sets <- list(
     "DCQN,RMA,A+B,FLN" = list(col=niceCols[6], name="QN,RMA,PM-only"),
     "DCQN,MIN1(PM-MM),RMA,A+B,FLN" = list(col=niceCols[2], name="QN,RMA,MIN(PM-MM,1)")
    );
  } else if (case == "RBC") {
    sets <- list(
     "QN,RMA,A+B" = list(col=niceCols[6], name="QN,RMA"),
     "RBC,QN,RMA,A+B" = list(col=niceCols[2], name="RBC,QN,RMA")
    );
  } else if (case == "ACC") {
    sets <- list(
     "QN,RMA,A+B,FLN,-X" = list(col=niceCols[2], name="QN,RMA,FLN"),
     "RMA,A+B,FLN,-X" = list(col=niceCols[2], name="RMA,FLN"),
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[6], name="ACC,RMA,FLN")
    );
  } else if (case == "RMAvMBEI") {
    sets <- list(
     "QN,RMA,A+B,FLN,-X" = list(col=niceCols[6], name="QN,RMA,FLN,-X"),
     "QN,MBEI,A+B,FLN,-X" = list(col=niceCols[2], name="QN,MBEI,FLN,-X"),
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[6], name="ACC,+300,RMA,FLN,-X"),
     "ACC,+300,MBEI,A+B,FLN,-X" = list(col=niceCols[2], name="ACC,+300,MBEI,FLN,-X")
    );
  } else if (case == 2) {
  } else {
    stop("Unknown case");
  }
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

# The indices of the reference set
if (refSet == "all") {
  references <- seq_len(ncol(truth));
} else if (refSet == "females") {
  references <- females;
} else if (refSet == "males") {
  references <- males;
}

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
    filename <- sprintf("%s.eps", imgName);
    eps(filename, width=5, height=5);
    par(mar=c(4,5,1,1)+0.1);
    par(cex=1, mgp=c(2.7,0.7,0));
  } else if (imgFormat == "screen") {
    par(mar=c(5,5,2,1)+0.1);
    par(cex=1, mgp=c(2.5,0.2,0));
  }
  par(cex.axis=1, cex.lab=1, font.lab=2);
  par(xaxs="i", yaxs="i");
  plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, las=1);
  abline(h=seq(0,1,by=by), lty=3, lwd=1);
  abline(v=seq(0,1,by=by), lty=3, lwd=1);
#  legend("bottomright", pch=19, lty=ltys, lwd=4, col=cols, legend=labels, bg="#eeeeee");
  for (kk in seq(along=sets)) {
    print(labels[kk]);
  
    set <- sets[[kk]];
    roc <- set$roc;
    plot(roc, line=TRUE, lwd=3, lty=ltys[kk], col=cols[kk], add=TRUE);
  } # for (kk in ...)
#  legend("bottomright", pch=19, lty=ltys, lwd=4, col=cols, legend=labels, bg="#eeeeee");
  legend("topleft", pch=19, lty=ltys, lwd=4, col=cols, legend=labels, bg="#eeeeee");

  if (imgFormat != "screen") {
    dev.off();
  }
}


if (case != "optimal")
  stop();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Graphs illustrating distribution of TP rate of single SNPs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (TRUE && !Device$isOpen(fig <- fig + 1)) {
  Device$set(fig);

  Device$setStyle("PowerPoint");
  par(mar=c(5,5,2,1)+0.1);
  par(cex=1, cex.axis=1, cex.lab=1,font.lab=2, mgp=c(2.5,0.2,0));
}



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
# Graphs illustrating ROC for a single SNP
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
name <- "ACC,+300,RMA,A+B,FLN,-X";
doROCone <- (name %in% names(sets));
if (doROCone) {
}

jj <- 1000; # Normal
jj <- 1002; # Perfect separation
jj <- 107; # Poor separation
jj <- 118; # Really bad separation

if (doROCone && !Device$isOpen(fig <- fig + 1)) {
  Device$set(fig, height=4.5);

  # Get data for one SNP
  M <- sets2[[name]]$M;
  Mj <- M[jj,];
  tj <- 1-unname(truth[jj,]);
  n <- length(Mj);

  # Calls
  for (threshold in c(NA, -1.2, -0.25, -0.15, 0, 0.25)) {
    if (!is.na(threshold)) {
      cj <- c(0,1)[as.integer(Mj < threshold)+1];
      tp <- sum(cj == 1 & tj == 1);
      tn <- sum(cj == 0 & tj == 0);
      fp <- sum(cj == 1 & tj == 0);
      fn <- sum(cj == 0 & tj == 1);
      nmales <- sum(tj == 1);
      nfemales <- sum(tj == 0);
      sens <- tp / (tp+fn);  # = tp/#males
      spec <- 1- fp / (tn+fp);  # = fn/#females
      cat(sprintf("goodness=%.2f, badness=%.2f\n", sens, 1-spec));
    }
    
    Device$setStyle("PowerPoint");
    par(ask=TRUE);
    par(mar=c(5,5,3,1)+0.1);
    par(cex=1, cex.axis=1, cex.lab=1,font.lab=2, mgp=c(2.5,0.2,0));
    
    plotDensity(Mj[males], col="blue", xlab="M", xlim=c(-1.5,1.5), ylim=c(0,4));
    plotDensity(Mj[females], col="red", add=TRUE);
#    title(main=getUnitNames(cdf)[jj]);
  
    ylim <- par("usr")[3:4];
    # Truth
    y <- ylim[2] - 0.1*diff(ylim);
    ys <- rep(y,n);
    dy <- -0.02*diff(ylim);
    ys[females] <- ys[females] + dy;
    points(Mj, ys, pch=19, col=c("red", "blue")[tj+1]);
    text(mean(xlim), y, "truth", adj=c(0.5,-0.5));
  
    if (!is.na(threshold)) {
      # Calls
      y <- ylim[2] - 0.2*diff(ylim);
      ys <- rep(y,n);
      ys[females] <- ys[females] + dy;
      points(Mj, ys, pch=21, col=c("red", "blue")[cj+1]);
      text(mean(xlim), y, "calls", adj=c(0.5,-0.5));
    
      # Threshhold
      y <- 0.3*diff(ylim);
      text(threshold, y, "threshold", srt=90, adj=c(0.5,-0.5));
      abline(v=threshold, lty=2);
    
      legends <- c(
        sprintf("#TP=%d", tp),
        sprintf("#FP=%d", fp),
        sprintf("#TN=%d", tn),
        sprintf("#FN=%d", fn),
        "",
        sprintf("goodness=%.2f", sens),
        sprintf("badness=%.2f",1-spec)
      );
      legend("topright", legend=legends, bg="#eeeeee");  
    }
  }
}



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Density of TP rates across SNPs at fixed FP rate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!Device$isOpen(fig <- fig + 1)) {
  Device$set(fig, height=4.5);
  C <- truth[,cc];
  nbrOfFemales <- sum(truth[1,cc]);
  fpStep <- 1/nbrOfFemales;
#  fpRate <- 0.01;
  fpRate <- 1*fpStep;
#  fpRate <- 2*fpStep;
#  fpRate <- 3*fpStep;

  cols <- sapply(sets, FUN=function(set) set$col);
 
  for (fpRate in c(1,2,3)*fpStep/2) {
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
      filename <- sprintf("%s.%s", imgName, imgFormat);
      eps(filename, width=6, height=0.618*6);
      par(mar=c(4,5,1,1)+0.1);
      par(cex=0.7, mgp=c(2.7,0.7,0));
    } else if (imgFormat == "screen") {
      par(mar=c(5,5,2,1)+0.1);
      par(cex=1, mgp=c(2.5,0.2,0));
    }
    
    
    xlabTp <- "TP rates (%)";
    plotDensity(100*D, from=0, to=100, xlab=xlabTp, lwd=2, col=cols);
    stext(side=4, pos=0, sprintf("FP rates %.2f%%", 100*fpRate), cex=0.7);
    legend("topleft", legend=colnames(D), col=cols, pch=19);
    
    if (imgFormat != "screen") {
      dev.off();
    }

stop()
next;    
    
    if (imgFormat == "eps") {
      imgName <- sprintf("ROC,%s,TPcdf,fpRate%.4f", case, fpRate);
      imgName <- gsub("[.]", "_", imgName);
      filename <- sprintf("%s.%s", imgName, imgFormat);
      eps(filename, width=6, height=0.618*6);
      par(mar=c(4,5,1,1)+0.1);
      par(cex=0.7, mgp=c(2.7,0.7,0));
    } else if (imgFormat == "screen") {
      par(mar=c(5,5,2,1)+0.1);
      par(cex=1, mgp=c(2.5,0.2,0));
    }

    plotCumulativeDensity(100*D, spar=0.5, xlab=xlabTp, lwd=2, col=cols);
    stext(side=4, pos=0, sprintf("FP rates %.2f%%", 100*fpRate), cex=0.7);
    legend("topleft", legend=colnames(D), col=cols, pch=19);
    
    if (imgFormat != "screen") {
      dev.off();
    }


    if (imgFormat == "eps") {
      imgName <- sprintf("ROC,%s,TPcdfInv,fpRate%.4f", case, fpRate);
      imgName <- gsub("[.]", "_", imgName);
      filename <- sprintf("%s.%s", imgName, imgFormat);
      eps(filename, width=6, height=0.618*6);
      par(mar=c(4,5,1,1)+0.1);
      par(cex=0.7, mgp=c(2.7,0.7,0));
    } else if (imgFormat == "screen") {
      par(mar=c(5,5,2,1)+0.1);
      par(cex=1, mgp=c(2.5,0.2,0));
    }
    
    plotCumulativeDensity(100*(1-D), spar=0.5, xlab=xlabTp, lwd=2, col=cols, axes=FALSE);
    box();
    axis(side=2);
    at <- seq(0,100,by=20);
    axis(side=1, at=at, label=100-at);
    stext(side=4, pos=0, sprintf("FP rates %.2f%%", 100*fpRate), cex=0.7);
    legend("topleft", legend=colnames(D), col=cols, pch=19);
    
    if (imgFormat != "screen") {
      dev.off();
    }
  } # for (fpRate ...)
}

# Average performance of SNPs
Dm <- rowMedians.matrix(D, na.rm=TRUE);
#Dm <- apply(D, MARGIN=1, FUN=min, na.rm=TRUE);

# A SNP with 0% TP rate
jj000 <- which(Dm == 0);

# A SNP with 60% TP rate
jj600 <- which(0.59 <= Dm & Dm <= 0.61);

# A SNP with 80% TP rate
jj800 <- which(0.79 <= Dm & Dm <= 0.81);

# A SNP with 95% TP rate
jj999 <- which(Dm >= 0.999);


jj <- jj800;
jj <- setdiff(jj, 1:100)[1];

# Example: Moderate SNP
jj <- indexOf(cdf, "SNP_A-1790385");
unitNames <- getUnitNames(cdf);
jj <- which(unitNames == "SNP_A-1884486");  # POOR
jj <- which(unitNames == "SNP_A-1790385");  # MODERATE
snpName <- unitNames[jj];
print(c(jj, snpName));

if (!Device$isOpen(fig <- fig + 1)) {
  Device$set(fig);

  Device$setStyle("PowerPoint");

  if (imgFormat == "eps") {
    imgName <- sprintf("ROC,%s,%s", case, snpName);
    filename <- sprintf("%s.%s", imgName, imgFormat);
    eps(filename, width=6, height=6);
    par(mar=c(4,5,1,1)+0.1);
    par(cex=1, mgp=c(2.7,0.7,0));
  } else if (imgFormat == "screen") {
    par(mar=c(5,5,2,1)+0.1);
    par(cex=1, mgp=c(2.5,0.2,0));
  }

  par(cex.axis=1, cex.lab=1, font.lab=2);
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
  legend("bottomright", pch=19, lty=ltys, lwd=4, col=cols, legend=labels, bg="#eeeeee");

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
    plot(roc, line=TRUE, lwd=4, col=cols[kk], add=TRUE);

#    roc <- fitRoc(tj, Mj, call=toCall, ncuts=11);
#    plot(roc, pch=19, cex=1.2, show.thresh=TRUE, add=TRUE);
  }
  box();

  if (imgFormat != "screen") {
    dev.off();
  }
}

stop();



kk <- 2; 
s <- 1;
Mf <- sets[[kk]]$M[,females];
Mm <- sets[[kk]]$M[,males];
main <- names(r)[kk];

# Density plots males/females
dev.set(3);
par(mar=c(5,4,2,2)+0.1); 
plot(NA, xlim=c(-1,1)*2, ylim=c(0,2.5), xlab="M", ylab="Density", main=main);
lines(density(na.omit(Mm/s)), lwd=2, col="orange");
lines(density(na.omit(Mf/s)), lwd=2, col="purple"); 

# (x,M) plots males/females
dev.set(4);
par(mar=c(5,4,2,2)+0.1);
plot(NA, xlim=c(0,max(x)), ylim=c(-1,1)*2, ylab="M", xlab="x", main=main);
points(x,Mm[,1]/s, pch=19, cex=0.7, col="orange");
points(x,Mf[,1]/s, pch=19, cex=0.7, col="purple"); 


stop()

## # Add CNAT v4 results
## name <- "CNAT,QN,PLIER,FLN+GC,xport";
## if (!name %in% names(sets)) {
##   filename <- sprintf("%s,M23.RData", name);
##   pathname <- filePath(dataPath, filename);
##   M <- loadObject(pathname);
## 
##   # Reorder samples according to 'truth'.
##   arrays <- gsub("^exp_([0-9]*).cn$", "\\1", colnames(M));
##   arrays <- as.integer(arrays);
##   o <- order(arrays);
##   M <- M[,o];
## 
##   set <- list();
##   set$M <- M;
## 
##   sets[[name]] <- set;
## }

# Add D & V statistics too
cols <- c(D="black", V="blue");
srcName <- "ACC,DCQN,pm,-X,+300,RMA,A+B,FLN,-X";
srcName <- "ACC,+300,RMA,A+B,FLN,-X";
srcName <- "ACC,+350,RMA,A+B,FLN,-X";
for (what in c("D", "V")) {
  name <- paste(srcName, what, sep=",");
  if (!name %in% names(sets)) {
    set <- sets[[srcName]];
    s2 <- rowVars(set$M[,females], na.rm=TRUE);
    if (what == "D") {
      w <- 1/sqrt(s2);
    } else if (what == "V") {
      w <- 1/s2;
    }
    w <- w / median(w, na.rm=TRUE);
    set$M <- w * set$M;
    set$name <- paste(set$name, what, sep=",");
    set$theta <- NULL;
    set$roc <- NULL;
    set$auc <- NULL;
    set$col <- cols[what];
    sets[[name]] <- set;
  }
} # for (what ...)

srcName <- "ACC,+300,RMA,A+B,FLN,-X";
for (what in c("QN")) {
  name <- paste(srcName, what, sep=",");
  set <- sets[[srcName]];
  theta <- set$theta;
  theta[,males] <- normalizeQuantile(theta[,males]);  
  theta[,females] <- normalizeQuantile(theta[,females]);
  set$theta <- theta;
  set$M <- NULL;
  set$A <- NULL;
  set$roc <- NULL;
  set$auc <- NULL;
  set$col <- "black";
  sets[[name]] <- set;
}
