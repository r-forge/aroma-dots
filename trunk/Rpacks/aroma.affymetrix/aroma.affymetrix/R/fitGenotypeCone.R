setMethodS3("fitGenotypeCone", "matrix", function(y, alpha=c(0.10, 0.075, 0.05, 0.03, 0.01), q=2, Q=98, ...) {
  require(sfit) || throw("Package 'sfit' not found.");

  # Fit simplex of (y_A,y_B)
  fit <- cfit(y, alpha=alpha, q=q, Q=Q, ...);
  M <- fit$M;
  colnames(M) <- c("A", "B");
  clazz <- class(M);

  # Re-arrange vertices in the order (origin, AA, BB)
  M <- fit$M
  idxOrigin <- which.min(apply(M, MARGIN=1, FUN=function(u) sum(u^2)));
  origin <- M[idxOrigin,];
  M <- M[-idxOrigin,];
  idxBBAA <- order(apply(M, MARGIN=1, FUN=function(u) diff(u)));
  M <- M[idxBBAA,];
  M <- rbind(origin, M);
  rownames(M) <- c("origin", "AA", "BB");
  class(M) <- clazz;
  fit$M <- M;

  W <- M[c("AA","BB"),];
  W <- t(W) - origin;
  W <- W / W[1,1];
 
  # Find the inverse
  Winv <- solve(W);

  Minv <- t(Winv %*% (t(M)-origin));
  class(Minv) <- clazz;
  fit$Minv <- Minv;

  W <- t(W);
  Winv <- t(Winv);

  fit$origin <- origin;
  fit$W <- W;
  fit$Winv <- Winv;

  # Re-arrange X too
  if (!is.null(fit$X)) {
    fit$X <- fit$X[,oB];
  }

  fit$params <- list(
    alpha=alpha,
    q=q,
    Q=Q,
    ...
  )

  fit;
}) # fitGenotypeCone()


############################################################################
# HISTORY:
# 2006-05-08
# o Created.
############################################################################
