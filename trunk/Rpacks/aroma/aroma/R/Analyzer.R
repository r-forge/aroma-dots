setConstructorS3("DesignMatrix", function(matrix=NULL) {
  if (!is.null(matrix))
    matrix <- as.matrix(matrix);

  extend(Object(), "DesignMatrix",
    .matrix = matrix
  )
})

setMethodS3("as.matrix", "DesignMatrix", function(x) {
  # To please R CMD check...
  this <- x;

  this$.matrix;
})



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Analyzer
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
setConstructorS3("Analyzer", function(mao=NULL, design=NULL) {
  if (!is.null(mao)) {
    if (!inherits(mao, "MicroarrayData"))
      throw("Object to be analyzed is not a MicroarrayData object: ", data.class(mao));
    if (is.null(design))
      design <- matrix(1, nrow=nbrOfSlides(mao), ncol=1);
  }

  extend(Object(), "Analyzer", 
    .mao    = mao,
    .design = design
  )
})


setMethodS3("getDesignMatrix", "Analyzer", function(this) {
  this$.design;
})

setMethodS3("setDesignMatrix", "Analyzer", function(this, design) {
  this$.design <- design;
})



setMethodS3("qqplotTTest", "Analyzer", function(this, main="t-test diagnostic", col="black", pch=".", low=-5, high=+5) {
  tma <- this$tma;
  if (is.null(tma))
    ttest(this);
  opar <- par("mar");
  on.exit(par(opar));
  subplots(2, nrow=2);
  t <- unlist(tma$tstat[,])
  low  <- which(t <= low);
  high <- which(t >= high);
  hist(t, xlab="t", nclass=100, main="Histogram and quantile-quantile plot of t-statistics", col=9, cex=0.8);
  qqnorm(tma, "tstat", xlab="Quantiles of standard normal", ylab="t", pch=".")
  highlight(tma, low, col="red")
  highlight(tma, high, col="blue")
})



setMethodS3("plotTTest", "Analyzer", function(this, main="t-test diagnostic", col="black", pch=".", low=-5, high=+5) {
  if (is.null(this$tma))
    ttest(this);
  opar <- par("mar");
  on.exit(par(opar));
  subplots(4);
  tma <- this$tma;
  low  <- which(tma$tstat[,] <= low);
  high <- which(tma$tstat[,] >= high);
  tma$tnum <- tma$tstat * tma$stderr;
  tma$tnumAbs <- abs(tma$tnum);
  opar <- par(mar=c(4,5,4,2));
  plot(tma, "tstatvsAavg",     xlab=expression(bar(A)), ylab=expression(T==(bar(M[1])-bar(M[2]))/SE), col=col, pch=pch)
  highlight(tma, low, col="red")
  highlight(tma, high, col="blue")
  opar <- par(mar=c(4,5,4,2));
  plot(tma, "stderrvsAavg",    xlab=expression(bar(A)), ylab="SE", col=col, pch=pch)
  highlight(tma, low, col="red")
  highlight(tma, high, col="blue")
  opar <- par(mar=c(5,5,1,2));
  plot(tma, "tnumAbsvsAavg",   xlab=expression(bar(A)), ylab=expression(abs(bar(M[1])-bar(M[2]))), col=col, pch=pch)
  highlight(tma, low, col="red")
  highlight(tma, high, col="blue")
  opar <- par(mar=c(5,5,1,2));
  plot(tma, "stderrvstnumAbs", xlab=expression(abs(bar(M[1])-bar(M[2]))), ylab="SE", col=col, pch=pch)
  highlight(tma, low, col="red")
  highlight(tma, high, col="blue")
  subplots(1);
  mtext(main, line=-1, cex=par("cex.main"))
})


setMethodS3("ttest", "Analyzer", function(this, treatments=getTreatments(this$.mao), var.equal=TRUE, verbose=TRUE) {
  if (!is.null(treatments)) {
    if (length(treatments) != nbrOfSlides(this$.mao))
      throw("Argument 'treatments' must be of the same length as the number of slides in the MicroarrayData object: ", length(treatments));
    nbrOfTreatments <- length(unique(treatments));
    if (nbrOfTreatments != 1 && nbrOfTreatments != 2)
      throw("To do either a one-sample or two-sample t-test there can only be one or two different treatments, respectively: ", nbrOfTreatments);
    uniqueTreatments <- unique(treatments);
    class <- list(
      which(treatments == uniqueTreatments[1]),
      which(treatments == uniqueTreatments[2])
    );
  } else {
    nbrOfTreatments <- 1;
  }

  layout <- getLayout(this$.mao);
  genes <- getGeneGroups(layout);
  geneSpots <- getSpots(genes);
  nbrOfGenes <- length(geneSpots);
  nbrOfSpots <- nbrOfSpots(layout);

  # Calculate the t-statistic *genewise* and set the results to each spot.
  # Unfortunately this will be a bit slow, but it is the only way.
  Mavg <- matrix(NA, nrow=nbrOfGenes, ncol=nbrOfTreatments);
  Aavg <- matrix(NA, nrow=nbrOfGenes, ncol=nbrOfTreatments);
  vx   <- matrix(NA, nrow=nbrOfGenes, ncol=nbrOfTreatments);
  nx   <- matrix(NA, nrow=nbrOfGenes, ncol=nbrOfTreatments);

  if (nbrOfTreatments == 1) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # One-sample t-test
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (verbose)
      cat("Performing a one-sample t-test...");

    # Copy the data to local variables
    M <- this$.mao$M;
    A <- this$.mao$A;

    # Do the t-test for each gene
    for (gene in seq(geneSpots)) {
      # 0. Get the spots belonging to the gene
      spot <- geneSpots[[gene]];

      Mi <- as.vector(M[spot,]);
      Ai <- as.vector(A[spot,]);
      ok <- (!is.na(Mi) & !is.na(Ai));
      Mi <- Mi[ok];
      Ai <- Ai[ok];

      # Step 1 - single-sample t-statistics
      nx[gene,] <- length(Mi);
      Mavg[gene,] <- mean(Mi);
      vx[gene,] <- var(Mi);
  
      # Step 2 - mean A values
      Aavg[gene,] <- mean(Ai);
    }

    # Calculate the T and the SE columnwise.
    df <- nx-1;
    stderr <- sqrt(vx/nx);
    tstat <- Mavg/stderr;    #  T = mean(X)/SE(X) = mean(X)/(var(X)/n)

    if (verbose)
      cat("ok\n");
  } else {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Two-sample t-test
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (verbose)
      cat("Performing a two-sample t-test...");

    # Copy the data to local variables
    M1 <- this$.mao$M[,class[[1]]];
    A1 <- this$.mao$A[,class[[1]]];
    M2 <- this$.mao$M[,class[[2]]];
    A2 <- this$.mao$A[,class[[2]]];

    # Do the t-test for each gene and for each treatment group
    for (gene in seq(geneSpots)) {
      # 0. Get the spots belonging to the gene
      spot <- geneSpots[[gene]];

      # 1a. For the data from class 1, keep only non-NA values
      Mi1 <- as.vector(M1[spot,]);
      Ai1 <- as.vector(A1[spot,]);
      idx <- (!is.na(Mi1) & !is.na(Ai1));
      Mi1 <- Mi1[idx];
      Ai1 <- Ai1[idx];

      # 1b. For the data from class 1, keep only non-NA values
      Mi2 <- as.vector(M2[spot,]);
      Ai2 <- as.vector(A2[spot,]);
      idx <- (!is.na(Mi2) & !is.na(Ai2));
      Mi2 <- Mi2[idx];
      Ai2 <- Ai2[idx];
   
      # 2. Two-sample t-statistics
      nx[gene,] <- c(length(Mi1), length(Mi2));
      Mavg[gene,] <- c(mean(Mi1)  , mean(Mi2)  );
      vx[gene,] <- c(var(Mi1)   , var(Mi2)   );
  
      # 3. Mean A values
      Aavg[gene,] <- c(mean(Ai1)  , mean(Ai2)  );
    } # for (gene ...)

    # 4. Calculate the degrees of freedom and the standard error
    #    for each gene.
    if (var.equal) {
      # Two sample t-test
      df <- nx[,1] + nx[,2] - 2;
      v <- ((nx[,1]-1) * vx[,1] + (nx[,2]-1) * vx[,2]) / df;
      stderr <- sqrt(v * (1/nx[,1] + 1/nx[,2]));
    } else { 
      # Welch
      stderr1 <- sqrt(vx[,1]/nx[,1]);
      stderr2 <- sqrt(vx[,2]/nx[,2]);
      stderr  <- sqrt(stderr1^2 + stderr2^2);
      df      <- stderr^4/(stderr1^4/(nx[,1]-1) + stderr2^4/(nx[,2]-1));
    }

    # 5. Finally, calculate the t statistics for each gene.
    tstat <- (Mavg[,1] - Mavg[,2]) / stderr;

    if (verbose)
      cat("ok\n");
  } # if (nbrOfTreatments == ...)

  tma        <- MicroarrayData(layout=layout);
  tma$tstat  <- GeneSlideArray(tstat,  layout=layout);
  tma$df     <- GeneSlideArray(df,     layout=layout);
  tma$stderr <- GeneSlideArray(stderr, layout=layout);
  tma$Mavg   <- GeneSlideArray(Mavg,   layout=layout);
  tma$Aavg   <- GeneSlideArray(Aavg,   layout=layout);
  tma$.fieldNames <- c("tstat", "df", "stderr", "Mavg", "Aavg");

  an$tma <- tma;

  tma;
})  # ttest()


setMethodS3("lmGenewise", "Analyzer", function(this, field="M", weights=NULL) {
  nbrOfSlides <- nbrOfSlides(this$.mao);

  design <- getDesignMatrix(this);

  nbeta <- ncol(design);

  data <- this$.mao[[field]];

  layout <- getLayout(this$.mao);
  genes <- getGeneGroups(layout);

  geneSpots <- getSpots(genes);

  nbrOfGenes <- length(geneSpots);
  stdev.unscaled <- beta <- matrix(NA, nrow=nbrOfGenes, ncol=nbeta, 
                                     dimnames=list(NULL, colnames(design)));
  sigma <- rep(NA, length.out=nbrOfGenes);
  df.residual <- rep(0, length.out=nbrOfGenes);

  for (k in seq(along=geneSpots)) {
    idx <- geneSpots[[k]];
    y <- data[idx,];
    ok <- is.finite(y);
    y <- y[ok];
    idxLen <- length(idx);
    if (idxLen > 1 && nrow(design) != idxLen)
      des <- design %x% rep(1, idxLen)
    else
      des <- design;
    X <- des[ok,, drop=FALSE];
    if (!is.null(weights)) {
      w <- as.vector(weights[idx,][ok]);
    }

    if (length(y) >= nbeta) {
        if (is.null(weights))
          fit <- lm.fit(X, y)
        else
          fit <- lm.wfit(X, y, w);
      beta[k, ] <- fit$coef;

      stdev.unscaled[k, ] <- sqrt(diag(chol2inv(fit$qr$qr)));

      df.residual[k] <- fit$df.residual;
      if (df.residual[k] > 0) {
        if (is.null(weights))
          sigma[k] <- sqrt(sum(fit$residuals^2)/fit$df.residual)
        else
          sigma[k] <- sqrt(sum(w * fit$residuals^2)/fit$df.residual);
      }
    } # if (length(y) >= nbeta)
  } # for(...)
 
  fit <- list(coefficients=drop(beta), stdev.unscaled=drop(stdev.unscaled),
                                       sigma=sigma, df.residual=df.residual);
  invisible(fit);
})  # lmGenewise()



############################################################################
# HISTORY:
# 2002-11-20
# o Added lmGenewise().
# 2002-11-05
# o First successful use of GeneSlideArray. All the MicroarrayArray classes
#   can access either the [spot,slide] pair or the [gene,replicate,slide]
#   triplet, only based on the how many indices are used. This means that
#   any underlying plot functions etc do not have to worry. Yeah, 
#   polymorphism rules!
# 2002-10-29
# o Created!
############################################################################

