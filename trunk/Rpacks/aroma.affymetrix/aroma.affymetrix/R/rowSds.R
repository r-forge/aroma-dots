rowVars <- function(X, mean=NULL, ...) {
  n <- dim(X)[2];
  if (is.null(mean))
    mean <- rowMeans(X, ...);
  X <- (X - mean)^2;
  X <- rowMeans(X, ...);
  (n/(n-1)) * X;
}

rowSds <- function(X, ...) {
  sqrt(rowVars(X, ...));
}

