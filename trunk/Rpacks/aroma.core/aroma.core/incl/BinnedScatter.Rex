# Sample scatter data
n <- 10e3
x <- rnorm(n=n)
y <- rnorm(n=n)
xy <- cbind(x=x, y=sin(x)+y/5)

# Bin data and estimate densities
xyd <- binScatter(xy)

# Plot data
plot(xyd, pch=19)

# Thin scatter data by subsampling
rhos <- c(1/3, 1/4, 1/6)
for (kk in seq(along=rhos)) {
  xyd2 <- subsample(xyd, size=rhos[kk])
  points(xyd2, pch=19, col=kk+1)
}
