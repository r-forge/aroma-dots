#########################################################################/**
# @set "class=MAData"
# @RdocMethod setField
#
# @title "Sets one or many fields"
#
# @synopsis
#
# \arguments{
#   \item{value}{The values that should replace the specified subset of the
#     data. This argument should be a vector, matrix or a data frame and it
#     must have the same number of columns as the number of specified fields.}
#   \item{fields}{The field names to be set. If \@NULL, all fields are
#     set. Valid values are \code{"M"} and \code{"A"}.}
# }
#
# \description{
#   Sets one or many fields. The field names must be specified by the 
#   argument \code{field} or if \code{field} is \@NULL all fields are
#   considered. The argument \code{value} must be able to be converted to
#   a matrix by \code{as.matrix(value)} where the resulting matrix must
#   have the same number of columns as the specified number of fields.
# }
#
# \examples{
#   # The option 'dataset' is used to annotate plots.
#   options(dataset="sma:MouseArray")
#
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   # Set all log ratios (M) to zeros.
#   setField(ma, "M", rep(0, ma$size()))
#   range(extract(ma, M"))  # 0 0
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("setField", "MAData", function(this, fields, value) {
  value <- as.matrix(value);
  ncol <- nbrOfSlides(this);

  # df now contains the final settings. Now we just update field by field
  # without having to think about the arguments slides and index.
  k <- 1;
  for (field in fields) {
    if (field == "M") {
      this$M <- matrix(value[,k], ncol=ncol)
    } else if (field == "A") {
      this$A <- matrix(value[,k], ncol=ncol)
    } else {
      throw("Trying to set protected field. Fields that can be set are \"M\" and \"A\": ", field);
    }
    k <- k + 1;
  }
  invisible(this);
}, private=TRUE, deprecated=TRUE);


############################################################################
# TO BE REMOVED?!?
############################################################################
setMethodS3("getSlide", "MAData", function(this, slide) {
  if (slide < 1 || slide > nbrOfSlides(this))
    throw("Argument 'slide' is out of range: ", slide);
  M <- this$M[,slide];
  A <- this$A[,slide];
  slideName <- getSlideNames(this, slide=slide);
  colnames(M) <- slideName;
  colnames(A) <- slideName;
  MAData(M=M, A=A, layout=getLayout(this), extras=this$.extras)
}, trial=TRUE, deprecated=TRUE)

setMethodS3("read", "MAData", function(this, filename, path=NULL, layout=NULL, verbose=FALSE) {
  fields <- c("M", "A");
  res <- MicroarrayData$readToList(filename, path=path,
                                   reqFields=fields, verbose=verbose);
  
  # Create the MAData object.
  MAData(M=res$M, A=res$A, layout=layout)
}, static=TRUE, trial=TRUE);


############################################################################
# O B S O L E T E
############################################################################
setMethodS3("getGenewiseResiduals", "MAData", function(this, robust=TRUE, na.rm=TRUE) {
  res <- clone(this);

  fields <- c("M", "A");
  for (field in fields) {
    x <- this[[field]];
    z <- apply(x, MARGIN=1, FUN=zscore, robust=robust, na.rm=na.rm);
    z <- t(z);
    res[[field]] <- z;
  }

  setSlideNames(res, names=getSlideNames(this));
  res;
}, deprecated=TRUE)



#########################################################################/**
# @RdocMethod getGeneResiduals
#
# @title "Gets the genewise residuals"
#
# @synopsis
#
# \arguments{
#   \item{robust}{If \@TRUE the median will be used for calculating
#     the residuals, otherwise the mean will be used.}
# }
#
# \description{
#   Calculates the genewise residuals. What spots belongs to which genes is
#   defined by the layout of the slides, where all slides are assumed to
#   have the same layout. See class @see "Layout" for more 
#   information about this.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#   ma.norm <- clone(ma)
#   normalizeWithinSlide(ma.norm, method="s")
#   normalizeAcrossSlides(ma.norm)
#   res <- getGeneResiduals(ma)
#   res.norm <- getGeneResiduals(ma.norm)
#
#   # Statistics
#   print(summary(as.vector(res)))
#   print(summary(as.vector(res.norm)))  # Improvement
#   nimproved <- sum(abs(res) >= abs(res.norm), na.rm=TRUE)
#   cat("Number of 'improved' spots:", nimproved, "\n")
#   cat("Number of 'worsened' spots:", length(res)-nimproved, "\n")
#
#   # Plots
#   xlim <- range(res, res.norm, na.rm=TRUE);
#   subplots(4)
#   plot(ma); hist(res, nclass=50, xlim=xlim);
#   plot(ma.norm); hist(res.norm, nclass=50, xlim=xlim);
# }
#
# @author
#
# \seealso{
#   @seemethod "getGeneVariability"
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getGeneResiduals", "MAData", function(this, robust=TRUE) {
  layout <- getLayout(this);
  geneGroups <- getGeneGroups(layout);
  genes <- getSpots(geneGroups);
  ngenes <- nbrOfGroups(geneGroups);
  
  data  <- this$M;
  
  variability <- rep(0, ngenes);
  residuals <- data;
  if (robust == TRUE) {
    # Robust standardization
    for (l in 1:ngenes) {
      spots <- genes[[l]];
      x <- data[spots,];
      ok <- !is.na(x);
      x <- x - median(x[ok]);
      variability[l] <- 1.4826*median(abs(x[ok]));
      residuals[spots,] <- x;
    }
  } else {
    # Non-robust standardization
    for (l in 1:ngenes) {
      spots <- genes[[l]];
      x <- data[spots,];
      ok <- !is.na(x);
      x <- x - mean(x[ok]);
      variability[l] <- sqrt(var(abs(x[ok]))); # Standard Deviation
      residuals[spots,] <- x;
    }
  } # if (robust == ...)

  names(variability) <- names(genes);
  attr(variability, "df") <- as.numeric(getSizes(geneGroups));
  class(variability) <- "GeneVariability";
  attr(residuals, "discrepancies") <- variability;
  class(residuals) <- "GeneResiduals";
  residuals;
}, private=TRUE, deprecated=TRUE)


setMethodS3("getGeneDiscrepancies", "MAData", function(this, robust=TRUE, standardize=FALSE) {
  warning("getGeneDiscrepancies() is deprecated. Use getGeneVariability() instead.")
  getGeneVariability(this, robust=robust);
}, deprecated=TRUE)


#########################################################################/**
# @RdocMethod as.BMAData
#
# @title "Calculates the genewise b-statistics"
#
# \description{
#  Calculates the genewise b-statistics (B) according to [1].
#  This method relies on the sma function @see "sma::stat.bayesian",
#  which was implemented by the authors of [1]. However, this method might
#  be somewhat faster.
#
#  \bold{Note:} We recommend that you use Gordon Smyth's limma package 
#  (\url{http://bioinf.wehi.edu.au/limma/}) and the @see "limma::ebayes" 
#  function instead, which is also based on [1]. Moreover, the limma package 
#  is richly documented. /July, 2003.
# }
#
# @synopsis
#
# \value{Returns a @see "BMAData" object.}
#
# \examples{\dontrun{
#   Use ebayes() in G. Smyth's limma package instead.
# }}
#
# @author
#
# \references{
#  [1] Lönnstedt, I. and Speed, T. P. (2002). Replicated microarray data. 
#      \emph{Statistica Sinica} \bold{12}, 31-46. 
# }
#
# \seealso{
#   @see "sma::stat.bayesian".
#   @see "limma::ebayes" (\url{http://bioinf.wehi.edu.au/limma/}).
#   @see "BMAData".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("as.BMAData", "MAData", function(this, p=0.01, v=NULL, a=NULL, c=NULL, k=NULL, df=NULL, nw=1, useCache=TRUE, saveToCache=TRUE) {
  if (nbrOfSlides(this) < 2)
    throw("To calculate the b-statistics there need to be at least two slides.");

  params <- list(p=p, v=v, a=a, c=c, k=k);
  # if 'cache' is specified then 'X=this' is obsolete.
  if (useCache == TRUE && !is.null(B.Xprep <- this$.cache$B.Xprep)) {
    bayesian <- stat.bayesian(nb=df, nw=nw, Xprep=B.Xprep, para=params);
    rm(B.Xprep);
  } else {
    bayesian <- stat.bayesian(unclass(this$M), nb=df, nw=nw, para=params);
    if (saveToCache) {
      this$.cache$B.Xprep <- bayesian$Xprep;
    }
  }
  
  # These two are redundant, but useful until implemented my own version.
  #  cache$Mbar <- NULL;
  #  cache$nb <- NULL;
  A.mean <- apply(unclass(this$A), MARGIN=1, FUN=mean, na.rm=TRUE);
  bma <- BMAData(M=bayesian$Xprep$Mbar, A=A.mean,
                     B=bayesian$lods, df=bayesian$Xprep$nb,
                     params = bayesian$para,
                     layout=getLayout(this), extras=this$.extras)
  bma;
}, deprecated=TRUE)  # as.BMAData()


#########################################################################/**
# @RdocMethod as.TMAData
#
# @title "Calculates the genewise t-statistics"
#
# \description{
#  Calculates the genewise t-statistics (T) and also the number of
#  \emph{degrees of freedom} (df) for each gene, the mean log-ratio (M) and
#  the mean log-intensity (A) for each gene.
#  The t-statistics is defined as 
#  \eqn{T=\frac{\overline{M}}{SE(M)}}{T=mean(M)/SE(M)}, where \eqn{SE(M)} is
#  the standard error of the mean.
#  Note that calculated values of T, df, M and A are stored in a \emph{n}x1
#  matrix, where \emph{n} is the number of \emph{spots}. This means that
#  if there are \emph{within-slide} replicates, the values are actually
#  repeatedly stored. 
# }
#
# @synopsis
#
# \value{Returns a @see "TMAData" object.}
#
# \examples{
#   # The option 'dataset' is used to annotate plots.
#   options(dataset="sma:MouseArray")
#
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw, bgSubtract=TRUE)
#   tma <- as.TMAData(ma)
#
#   genes <- getGeneGroups(layout)
#   spots <- getSpots(genes)
#   idx <- unlist(lapply(spots, FUN=function(x) x[1]))
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("as.TMAData", "MAData", function(this, treatments=getTreatments(this), var.equal=TRUE, verbose=TRUE) {
  if (!is.null(treatments)) {
    if (length(treatments) != nbrOfSlides(this))
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

  layout <- getLayout(this);
  genes <- getGeneGroups(layout);
  geneSpots <- getSpots(genes);
  nbrOfGenes <- length(geneSpots);
  nbrOfSpots <- nbrOfSpots(layout);

  # Calculate the t-statistic *genewise* and set the results to each spot.
  # Unfortunately this will be a bit slow, but it is the only way.
  Mx <- matrix(NA, nrow=nbrOfSpots, ncol=nbrOfTreatments);
  vx <- matrix(NA, nrow=nbrOfSpots, ncol=nbrOfTreatments);
  Ax <- matrix(NA, nrow=nbrOfSpots, ncol=nbrOfTreatments);
  nx <- matrix(NA, nrow=nbrOfSpots, ncol=nbrOfTreatments);

  if (nbrOfTreatments == 1) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # One-sample t-test
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (verbose)
      cat("Performing a one-sample t-test...");

    # Copy the data to local variables
    M <- unclass(this$M);
    A <- unclass(this$A);

    # Do the t-test for each gene
    for (spot in geneSpots) {
      Mi <- as.vector(M[spot,]);
      Ai <- as.vector(A[spot,]);
      ok <- (!is.na(Mi) & !is.na(Ai));
      Mi <- Mi[ok];
      Ai <- Ai[ok];

      # Step 1 - single-sample t-statistics
      nx[spot,] <- length(Mi);
      if (any(nx[spot,] > 1)) { # because if Mi is empty var() gives an error.
        Mx[spot,] <- mean(Mi);
        vx[spot,] <- var(Mi);
      }
  
      # Step 2 - mean A values
      Ax[spot,] <- mean(Ai);
    }

    # Calculate the T and the SE columnwise.
    df <- nx-1;
    stderr <- sqrt(vx/nx);
    tstat <- Mx/stderr;    #  T = mean(X)/SE(X) = mean(X)/(var(X)/n)

    if (verbose)
      cat("ok\n");
  } else {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Two-sample t-test
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (verbose)
      cat("Performing a two-sample t-test...");

    # Copy the data to local variables
    M1 <- unclass(this$M[,class[[1]]]);
    A1 <- unclass(this$A[,class[[1]]]);
    M2 <- unclass(this$M[,class[[2]]]);
    A2 <- unclass(this$A[,class[[2]]]);

    # Do the t-test for each gene and for each treatment group
    for (spot in geneSpots) {
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
      nx[spot,] <- c(length(Mi1), length(Mi2));
      if (any(nx[spot,1] > 1)) { # because if Mi1 is empty var() gives an error.
        Mx[spot,1] <- mean(Mi1);
        vx[spot,1] <- var(Mi1);
      }
      if (any(nx[spot,2] > 1)) { # because if Mi2 is empty var() gives an error.
        Mx[spot,2] <- mean(Mi2);
        vx[spot,2] <- var(Mi2);
      }
  
      # 3. Mean A values
      Ax[spot,] <- c(mean(Ai1)  , mean(Ai2)  );
    } # for (spot in geneSpots)

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
    tstat <- (Mx[,1] - Mx[,2]) / stderr;

    if (verbose)
      cat("ok\n");
  } # if (nbrOfTreatments == ...)

  TMAData(M=Mx, A=Ax, T=tstat, df=df, stderr=stderr, layout=layout, extras=this$.extras)
}, deprecated=TRUE)  # as.TMAData()


############################################################################
# HISTORY:
# 2006-02-08
# o Although deprecated, as.TMAData() still gave a warning about "...if 
#   (nx[spot, ] > 1) { : the condition has length > 1 and only the first 
#   element will be used".  Fixed.
# 2005-07-21
# o Made as.BMAData() and as.TMAData() deprecated.
# 2005-07-19
# o Replaced all path="" arguments to path=NULL.
# 2002-11-12
# o Created. Here all deprecated functions will be put in quarantine until
#   it is safe to remove them. 
############################################################################
