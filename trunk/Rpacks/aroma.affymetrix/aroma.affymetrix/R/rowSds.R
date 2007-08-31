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


rowMads <- function(X, centers=rowMedians(X, na.rm=na.rm), na.rm=FALSE, constant=1.4826, ...) {
  constant * rowMedians(abs(X - centers), na.rm=na.rm);
} # rowMads()
