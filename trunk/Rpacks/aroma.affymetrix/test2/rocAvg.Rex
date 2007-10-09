# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
savehistory()
library(aroma.affymetrix);
library(R.graphics);
library(R.native);  # rowMedians(..., na.rm=TRUE)
source("fitRoc.R");
log <- Arguments$getVerbose(-5);
timestampOn(log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
getBlockAverageMap <- function(n, h=1, ...) {
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
cc <- which(!colnames(truth) %in% excl)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

R <- 5;

case <- "optimal";
#case <- "ACC";
#case <- "steps";

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
     "CNAG,FLN+GC,r60,xport" = list(col=niceCols[8], name="CNAG*"),
     "CNAT,QN,PLIER,FLN+GC,xport" = list(col=niceCols[4], name="CNAT v4"),
     "dChip,ISN,MBEI,A+B,xport" = list(col=niceCols[2], name="dChip"),
#     "dChip,QN,MBEI,A+B,xport" = list(col=niceCols[1], name="dChip,QN"),
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[6], name="CRMA")
#     "ACC,+300,RMA,A+B,FLN,-X,GCN,-X" = list(col=niceCols[5], name="CRMA,GC"),
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
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[6], name="ACC,RMA,FLN")
    );
  } else if (case == "steps") {
    sets <- list(
     "QN,RMA,A+B" = list(col=niceCols[2], name="QN"),
     "ACC,+300,RMA,A+B" = list(col=niceCols[4], name="ACC"),
     "ACC,+300,RMA,A+B,FLN,-X" = list(col=niceCols[6], name="ACC+FLN")
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
  } else if (case == "allVdiploid") {
    sets <- list(
     "ACC,+300,RMA,A+B,FLN,-X;all" = list(col=niceCols[6], name="CRMA,all"),
     "ACC,+300,RMA,A+B,FLN,-X;adj" = list(col=niceCols[7], name="CRMA,adj"),
     "ACC,+300,RMA,A+B,FLN,-X;diploid" = list(col=niceCols[5], name="CRMA,diploid"),
     "dChip,ISN,MBEI,A+B,xport;adj" = list(col=niceCols[1], name="dChip,adj"),
     "dChip,ISN,MBEI,A+B,xport;diploid" = list(col=niceCols[2], name="dChip,diploid")
    );
  } else if (case == 2) {
  } else {
    stop("Unknown case");
  }
}



# The indices of the reference set
if (refSet == "all") {
  references <- seq_len(ncol(truth));
} else if (refSet == "females") {
  references <- females;
} else if (refSet == "males") {
  references <- males;
}


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
for (R in seq(from=1, to=10, by=0.5)) {
  if (R == 1)
    next;
print(R);

idxs <- getBlockAverageMap(n=nrow(truth), h=R);
Rapprox <- attr(idxs, "hApprox");

verbose && printf(verbose, "R=%.3f, Rapprox=%.3f\n", R, Rapprox);

if (!identical(attr(x, "isOrdered"), TRUE)) {
  o <- order(x);
  truth <- truth[o,];
  x <- x[o];
  attr(x, "order") <- o;
  attr(x, "isOrdered") <- TRUE;
}


for (kk in seq(along=sets)) {
  name <- names(sets)[kk];
  if (regexpr("R=[0-9]", name) != -1)
    next;

  newName <- sprintf("%s,R=%.3f", name, R);
  pathname <- filePath(dataPath, sprintf("%s.RData", newName));
  if (isFile(pathname))
    next;
  print(pathname);
#  if (newName %in% names(sets))
#    next;


  verbose && enter(verbose, sprintf("Set #%d ('%s') of %d", 
                                                   kk, name, length(sets)));
  set <- sets[[name]];

  # Get reordered data
  o <- attr(x, "order");
  set$x <- blockAvg(x[o], idxs);
  set$M <- blockAvg(set$M[o,], idxs);
  if (!is.null(set$A))
    set$A <- blockAvg(set$A[o,], idxs);
  set$truth <- blockAvg(truth[o,], idxs);
  set$roc <- NULL;
  set$auc <- NULL;
  set$theta <- NULL;
  set$name <- sprintf("%s,R=%.3f", set$name, R);

  cat("Saving...\n");
  saveObject(set, pathname);
  cat("Saving...ok\n");

#  sets[[newName]] <- set;

  verbose && exit(verbose);
} # for (kk in ...)

} # for (R in ...)