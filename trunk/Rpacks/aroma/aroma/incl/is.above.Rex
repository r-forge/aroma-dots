  is.above <- function(line, xy, na.replace=FALSE) {
    x.idx <- apply(as.matrix(xy$x), MARGIN=1, FUN=function(x) {
        x0 <- which(x >= line$x);
        if (length(x0) == 0) 1 else max(x0);
      }
    )

    xs <- matrix(c(line$x[x.idx], line$x[x.idx+1]), nrow=2, byrow=TRUE); 
    delta <- (xy$x-xs[1,])/(xs[2,]-xs[1,]);

    ys <- matrix(c(line$y[x.idx], line$y[x.idx+1]), nrow=2, byrow=TRUE);

    y.at.x <- (ys[2,]-ys[1,])*delta+ys[1,];
    y.dist <- xy$y-y.at.x;

    res <- (y.dist > 0);
    if (length(res) > 0 && !is.null(na.replace)) {
      res[is.na(res)] <- na.replace;
    }
    res;
  }
