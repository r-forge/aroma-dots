source("init.R");

imgFormat <- c("x11", "eps", "png")[3];

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Do boxplots for all log-ratios
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
Mlim <- c(-2,1.5);
Mlab <- expression(M==log[2](theta/theta[R]));
for (kk in seq(along=sets)) {
  set <- sets[[kk]];
  verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", kk, set$name, length(sets)));

  name <- set$name;
  name <- gsub("[*]", "+", name);
  imgName <- paste(c(name, "M", "boxplot"), collapse=",");

  figPathT <- filePath(figPath, "boxplots");
  if (!devIsOpen(imgName)) {
    width <- imgScales[imgFormat] * 1.2;
    height <- 3/4*0.618*width;
    if (imgFormat == "x11") {
      devOpen(imgFormat, width=width, height=height, name=imgName);
    } else {
      mkdirs(figPathT);
      imgPathname <- filePath(figPathT, sprintf("%s.%s", imgName, imgFormat));
      if (!force && isFile(imgPathname)) {
        verbose && cat(verbose, "Already exists: ", imgPathname);
        verbose && exit(verbose);
        next;
      }
      devOpen(imgFormat, imgPathname, width=width, height=height, name=imgName);
    }

    set <- addRocData(set, verbose=verbose);
    M <- getData(set$rocData, ordered=FALSE, complete=FALSE);
    M <- as.data.frame(M);
    gc <- gc();
    xlim <- c(1/2,ncol(M)+1/2);
    par(mar=c(5,4,1,1)+0.1);
    plot(NA, xlim=xlim, ylim=Mlim, las=2, xlab="", ylab=Mlab, axes=FALSE);
    abline(h=0, lty=3, col="#999999");
    cols <- colors[c(5,9)];
    boxplot(M, col=cols[n23[colnames(M)]], outline=FALSE, las=2, add=TRUE);
    stext(side=3, pos=0, cex=0.7, set$name);
    rm(M);
    gc <- gc();
    
    devDone();
  }

  verbose && exit(verbose);
}
