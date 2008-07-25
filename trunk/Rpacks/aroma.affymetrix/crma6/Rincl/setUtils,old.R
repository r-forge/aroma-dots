addLogRatios <- function(set, ...) {
  # Already added?
  if (!is.null(set$M))
    return(set);

  cat("Adding log-ratios: ", set$name, "\n", sep="");
  set$M <- loadLogRatios(set, ...);
  gc <- gc();

  set;
} # addLogRatios()


addRoc <- function(set, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Already calculated?
  if (!is.null(set$rocFit))
    return(set);

  verbose && enter(verbose, "Calculating ROC");
  verbose && cat(verbose, "Data set: ", set$name);

  # Get log-ratios
  M <- set$M;

  # Get true CNs
  C <- set$truth;
  if (is.null(C)) {
    verbose && enter(verbose, "Creating default truth");
    # FIX?!?
    C <- matrix(as.integer(1), nrow=nrow(M), ncol=ncol(M));
    C[,(n23[colnames(M)] == 2)] <- as.integer(0);
    verbose && exit(verbose);
  }

#  subset <- 1:10e3; M <- M[subset,]; C <- C[subset,]; rm(subset);

  # Vectorize
  M <- as.vector(M);
  C <- as.vector(C);

print(table(C));
    
  # Keep only finite data points
  ok <- (is.finite(M) & is.finite(C));
  n <- length(M);
  callRate <- sum(ok)/n;
  cat(sprintf("Keeping %d out of %d (%.4f%%) data points, i.e. excluding %d (%.4f%%) \n", sum(ok), n, 100*callRate, sum(!ok), 100*sum(!ok)/n));

  M <- M[ok];
  C <- C[ok];
  rm(ok);
  gc <- gc();
  
  verbose && enter(verbose, "Fitting ROC");
  verbose && cat(verbose, "Log-ratios (M):");
  verbose && str(verbose, M);
  verbose && cat(verbose, "Class (C):");
  verbose && str(verbose, C);

  # fit ROC curve
  rocData <- RocData(truth=C, data=M, hasNA=FALSE);
  print(rocData);
  fit <- fit(rocData, ncuts=1200, verbose=less(verbose,1));
  print(fit);
  rm(C, M);
  verbose && exit(verbose);
 
  rocFit <- list(nbrOfLoci=n, callRate=callRate, fit=fit, auc=auc(fit));
  verbose && cat(verbose, "ROC fit:");
  verbose && str(verbose, rocFit);
  
  set$rocFit <- rocFit;
  rm(rocFit);

  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  set;
} # addRoc()
